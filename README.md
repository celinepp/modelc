# modelc
Monte Carlo Methods for Derivative Pricing
This C++ project involved researching and implementing Monte Carlo simulation techniques for pricing financial derivatives. The core focus was on developing robust models for path generation and pricing methods using C++.

Overview
The key components developed include:

Classes for stochastic process modeling such as Black-Scholes, Ornstein-Uhlenbeck, and square root diffusion.
Path generation frameworks using Brownian bridges and correlated random variates.
Valuing European options through variance reduction techniques.
An extensible architecture for adding new pricing models and methods.
The project showcases research skills in quantitative finance, C++ programming, and computational methods for derivative pricing.

Stochastic Process Modeling
Classes were created to model common stochastic processes used in derivative pricing models, including geometric Brownian motion, Ornstein-Uhlenbeck, and square root diffusion. This provides a library of process models that can be reused across pricing applications.

Path Generation
Efficient and accurate path generation is critical for Monte Carlo pricing. Techniques implemented include Brownian bridge construction and generating correlated random variates for multidimensional processes.

Variance Reduction
Variance reduction techniques like antithetic variates, control variates, and dynamic replication were implemented to improve the accuracy and performance of the pricing methods.

Extensible Design
The project uses an extensible architecture so new pricing models, path generators, and variance reduction techniques can be easily added. This enables ongoing expansion as research continues.

Next Steps
Future work involves expanding the model classes, path generators, and pricing engines. Areas of research include quasi-Monte Carlo methods, additional stochastic processes, and GPU/parallel computing.

The project demonstrates core skills in quantitative research, C++ programming, and computational finance. The analysis and implementation of these Monte Carlo methods and models provides a solid foundation for further derivatives pricing work.
