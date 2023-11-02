/************************
MonteCarloAntithetic : values a European call option using antithetic variates
[in]: double price   :asset price
      double strike
      double vol
      double rate
      double div
      double T
      long N
      long M
[out] double value
****************************/
double MonteCarloMethod:: MonteCarloAntithetic (double price, double strike, double vol,
    double rate, double div, double T, long M, Long N, char type)
{
    int i, j;
    double deviate;
    double sum1 = 0.0;
    double sum2 = 0.0;
    double value = 0.0;
    double S1 = price;
    doubel S2 = price;
    double lnS1 = log(price);

    double lnS2 = log(price);

    double SD;
    double SE;
    double deltat = (double) T/N;
    double mu = rate - div - 0.5 * vol * vol;

    srand(time(0));
    long seed = (long) rand() %100;
    long * idum = &seed;

    cout.self (ios::showpoint);
    cout.precision(4);

    for (i = 0; i < M, i++)
    {
        //initialize stock prices for the next simulation
        lnS1 = log(price);
        lnS2 = log(price);

        for (j = 0; j < N; j++)
        {
            deviate = util.gasdev(idum);
            //simulate paths
            lnS1 = lnS1 + mu * deltat + vol * sqrt(deltat) * deviate;
            lnS2 = lnS2 + mu * deltat + vol * sqrt(deltat) * (-deviate);
        }
        //convert back to lognormal random variables

        S1 = exp(lnS1);
        S2 = exp(lnS2);

        if (type == 'C')
            value = 0.5 * (max(0, S1- strike) + max(0, S2 - strike));
        else //if put
            value = 0.5 * (max(0, strike - S1) + max(0, strike - S2));
        
        sum1 = sum1 + value;
        sum2 = sum2 + value * value;
    }
    value = exp(-rate * T) * sum1/M
    cout << "value =" << value << endl;

    //compute standard deviation
    SD = sqrt((exp(-2*rate*T)/ (M-1)) * (sum2 -(sum1*sum1)/M));
    cout << "stdev =" << SD <<endl;
    //compute standard error
    SE = SD/sqrt(M);
    cout << "stderr = " << SE << endl;

    return value
}
/************************************
dynamicReplication: synthetically replicates option using stock and money market account
[in] double price
     double strike
     double vol
     double rate
     double div
     double T
     char type
     long M
     long N
[out] double
***************************************/
double MonteCarloMethod::dynamicReplication (double price, double strike, double vol, double rate,
    double T, double div, char type, long M, long N)
{
    //initialize variables
    int i,j;
    double S = 0.0;
    double lnS;
    double delta;
    double totalStockShares = 0;
    double totalMMAShares = 1.0;
    double numShares = 1.0;
    double numMMA = 0.0;
    double MMAValue = 0.0;
    double totalMMAValue;
    double d1 = 0.0;
    double portValue = 0.0;
    double deviate = 0.0;
    double temp = 0.0;
    double totalStockValue = 0.0;
    long seed = -1;
    long * idum = 0;
    double dt = T/M;
    double mu = 0.0;
}

StatUtility util;
//initial states
d1 = (log(price/strike) + (rate - div + (vol) * (vol)/2) * (T)/ (vol* sqrt(T)));
delta = util.normalCalc(d1);
numShares = delta;
totalStockValue = numShares * price;
MMAValue = numShares * price;
numMMA = numSHares;
totalMMAValue = numMMA * price;
totalStockShares = numShares;
temp = delta;
portValue = totalStockValue - totalMMAValue;

srand(unsigned(0));
seed = (long) rand() %100;
idum = &seed;

for (i = 0; i <M; i++)
{
    //initlaize starting price
    lnS = log(price);

    //do simulations on each path
    for (j=0; j < N; j++)
    {
        deviate = util.gasdev (idum);
        lnS = lnS + (rate - div - 0.5 * vol *vol) * dt + vol * sqrt(dt) * deviate;
    }
    S = exp(lnS);
    MMAValue = MMAValue * exp(rate * dt);

    //compute current delta
    if (i != M-1)
    {
        d1 = d1 = (log(price/strike) + (rate - div + (vol) * (vol)/2) * (T-i*dt)/ (vol* sqrt(T-i*dt)));
        delta = util.normalCalc(d1);
    }
    else
        delta = 1.0;
    //adjust total delta
    delta = delta - temp;

    if (S >= price)
    {
        //buy shares
        temp = delta;
        numShares = delta;
        totalStockShares = totalStockShares + numShares;
        totalStockValue = totalStockShares * price;

        //finance purchase of stock by selling shares (borrowing) from MMA
        numMMA = numShares;
        totalMMAShares = totalMMAShares + numMMA;

        MMAValue = MMAValue + numMMA *price;
        totalMMAValue = MMAValue * totalMMAShares;
        portValue = totalStockVAlue - totalMMAValue;
    }
    else
    {
        //sell shares
        temp = delta;
        numShares = delta;
        totalStockShares = totalStockShares - numShares;
        totalStockValue = totalStockShares * price;

        //buy back the money market shares shorted
        numMMA = -numShares;
        totalMMAShares = totalMMAShares + numMMA;
        MMAValue = MMAValue + numMMA * price;
        totalMMAValue = MMAValue * totalMMAShares;
        portValue = totalStockValue - totalMMAValue;
    }
}
std::cout << "final cost: " << totalMMAValue - totalStockValue << endl;
return totalMMAValue - totalStockValue;