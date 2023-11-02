#include "handle.h"
typedef size_t Size;
typedef double Time;

namespace QuantLib
{
    namespace MonteCarlo
    {
        //General purpose Monte Carlo model for path samples
        /*Three components
            - S, sample accumlator
            - PG, path generator,
            - PP, path pricer
        The minimal interfaces for classes S, PG, and PP:
        class S
        {
            void add(VALUE_TYPE sample, double weight) const;
        };
        class PG
        {
            Sample<PATH_TYPE> next() const;
        };
        class PP:: unary_function<PATH_TYPE, VALUE_TYPE>
        {
            VALUE_TYPE operator() (PATH_TYPE &) const;
        };
        */
    template<class S, class PG, class PP>
    class MonteCarloModel
    {
        public:
            typedef typename PG::sample_type sample_type;
            typedef typename PP::result_type result_type;
            MonteCarloModel (const Handle<PG>& pathGenerator,
            const Handle<PP>& pathPricer,
            const S& sampleAccumulator,
            const Handle<PP>& cvPathpricer = Handle<PP>(),
            result_type cvOptionValue  = result_type());
            void addSamples(Size samples);
            const S& sampleAccumulator(void) const;
        private:
            Handle<PG> pathGenerator_;
            Handle<PP> pathPricer_;
            S sampleAccumulator_;
            Handle<PP> cvPathPricer_;
            result_type cvOptionValue_;
            bool isControlVariate_;
    };
    //inline definitions
    template<class S, class PG, class PP>
    inline void MonteCarloModel <S, PG, PP>::addSamples(Size samples)
    {
        for(Size j = 1; j <= samples; j++)
        {
            sample_type path = pathGenerator_->next();
            result_type price = (*pathPricer_)(path.value);
            if (isControlVariate_)
                price += cvOptionValue_-(*cvPathPricer_)(path.value);
            sampleAccumulator_.add(price, path.weight);
        }
    }
    template<class S, class PG, class PP>
    inline const S& MonteCarloModel <S, PG, PP>::sampleAccumulator() const {
        return sampleAccumulator_;
        }
    }
}

#include "ql/MonteCarlo/path.h"
#include "ql/RandomNumbers/randomarraygenerator.h"
#include "ql/diffusionprocess.h"

namespace QuantLib{
    namespace MonteCarlo{
        /******************
         Generates random paths from a random number generator
        **********************/
       template <class RNG>
       class PathGenerator {
        public:
            typedef Sample<Path> sample_type;
            //constructors
            PathGenerator (double drift,
                double variance,
                Time length,
                Size timeSteps,
                long seed = 0);
            //warning the initial time is assumed to be zero
            //and must not be included in the passed vector
            PathGenerator(double drift,
                double variance,
                const std::vector <Time>& times,
                long seed = 0);
            PathGenerator (double drift,
                double variance,
                const std::vector <Time>& times,
                long seed = 0);
            PathGenerator(const std::vector <double>& drift,
                const std::vector<double>& variance,
                const std::vector<Time>& times,
                long seed = 0)
            private:
                mutable Sample<Path> next_;
                Handle<RandomNumbers::RandomArrayGenerator<RNG> > generator_;                
       };
       template <class RNG>
       PathGenerator<RNG>::PathGenerator (double drift, double variance, Time length, Size timeSteps, long seed):
       next_(Path(timeSteps), 1.0)
       {
        QL_REQUIRE(timeSteps > 0, "PathGenerator: Time steps("+
        IntegerFormatter::toString(timeSteps) +") must be greater than zero");
        QL_REQUIRE(length > 0, "PathGenerator: length must be > 0");

        Time dt = length/timeSteps;

        for (Size i = 0; i<timeSteps; i++)
        {
            next_.value.times()[i] = (i+1) * dt;
        }
        next_.value.drift() = Array(timeSteps, drift * dt);

        QL_REQUIRE(variance >= 0.0, "PathGenerator: negative variance");
        generator_= Handle<RandomNumbers::RandomArrayGenerator<RNG> . (
            new RandomNumbers::RandomArrayGenerator <RNG>(
                Array(timeSteps, variance * dt), seed));
       }
       template <class RNG>
       PathGenerator <RNG> ::PathGenerator(double drift, double variance,
        const std::vector<Time>& times, long seed)
        : next_(Path(times.size()), 1.0)
        {
            QL_REQUIRE(variance >= 0.0, "PathGenerator: negative variance");
            QL_REQUIRE(time.size() >0, "PathGenerator: no times given");
            QL_REQUIRE(times[0] >= 0.0, "PathGenerator: first time("+
                DoubleFormatter::toString(times[0])+") must be non negative");
            Array variancePerTime(times.size());
            Time dt = times[0];
            next_.value.drift()[0] = drift * dt;
            variancePerTime[0] = variance * dt;
            for (Size i = 1; i < times.size(); i++)
            {
                //check current time is greater than previous time
                QL_REQUIRE(times[i] >= times[i-1],
                    "MultiPathGenerator: time("+ IntergerFormatter::toString(i-1)+") = "
                    + DoubleFormatter::toString(times[i1])
                    "is later than time ("+ IntergerFormatter::toString(i) +")=" +
                    DoubleFormatter::toString(times[i]));

                    dt = times[i] - times[i-1];
                    next_.value.drift()[i] = drift * dt;
                    variancePerTime[i] = variance * dt;
            }
            next_.value.times() = times;

            generator_= Handle<RandomNumbers::RandomArrayGenerator<RNG> >(
                new RandomNumbers::RandomArrayGenerator<RNG>(variancePerTime, seed));
        }
        template<class RNG> PathGenerator<RNG>::PathGenerator(
            const std::vector<double>& drift,
            const std::vector<double>& variance,
            const std::vector<Time>& times, long seed) : next_(Path(times.size()), 1.0)
            {
                //data validity check
                QL_REQUIRE(times.size() > 0, "PathGenerator: no times given");
                QL_REQUIRE(times[0] >= 0.0, "PathGenerator: first time(" +
                    DoubleFormatter::toString(times[0])+") must be non negative");
                QL_REQUIRE(variance.size() == times.size(),
                    "PathGenerator: mismatch between variance and time arrays");
                QL_REQUIRE(drift.size() == times.size(),
                    "PathGenerator: mismatch between drift and time arrays");
                
                Array variancePerTime(times.size());
                double dt = times[0];
                next_.value.drift()[0] = drift[0] * dt;

                QL_REQUIRE(variance[0] >= 0.0, "PathGenerator: negative variance");
                variancePerTime[0] = variance[0]*dt;
                
                for(Size i =1; i < times.size(); i++)
                {
                    QL_REQUIRE(variance[i] >= 0.0, "PathGenerator: negative variance");
                    variancePerTime [i] = variance[i] *dt;
                }
                next_.value.times() = times;

                generator_ = Handle<RandomNumbers::RandomArrayGenerator<RNG> >(
                    new RandomNumbers::RandomArrayGenerator<RNG>(variancePerTime,
                    seed));       
            }
            template <class RNG> inline const typename PathGenerator<RNG> :: sample_type&
            PathGenerator<RNG>::next() const
            {
                const Sample<Array>& sample = generator_->next();
                next_.weight = sample.weight;
                next_.value.diffusion() = sample.value;
                return next_.;
            }
    }
}

namespace QuantLib
{
    namespace MonteCarlo
    {
        /************************
        Path class for handling computations of drift and diffusion terms along a path single factor random walk
        *************************/
       class Path
       {
        public:
        Path(Size size);
        Path(const std::vector<Time>& times, const Array&drift, const Arrays diffusion);
        //inspectors
        double operator[] (int i) const;
        Size size() const;
        //read/write access to components
        const std::vector<Time>& times() const;
        std::vector<Time>& times();
        const Array& drift() const;
        Array& drift();
        const Array& drift() const;
        Array& diffusion();
        private:
        std::vector<Time> times_; //vector of time instances
        Array drift_;
        Array diffusion_;
       };
       //inline definitions
       inline Path::Path (Size size)
        :times_(size), drift_(size), diffusion_(size) {}
        inline Path::Path(const std::vector<Time>& times, const Array& drift, const
            Array& diffusion)
            :time_(times), drift_(drift), diffusion_(diffusion)
        {
            QL_REQUIRE(drift_.size() == diffusion_.size();
                "Path: drift and diffusion have different size");
            QL_REQUIRE(times_.size() == drift_.size(),
                "Path: times and drift have different size");
        }
        //overload [] operator
        inline double Path::operator[] (int i) const
        {
            return drift_[i] + diffusion_[i];
        }
        inline Size Path::size() const 
    }
}