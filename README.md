### DTQbp: Modeling Detection, Contact Tracing and Quarantine with Branching Processes


We implemented a computational Branching Processes model for disease transmission, under control measures including case detection, contact tracing and quarantine. The main functions of this model are

<ol>
  <li>Projecting the epidemic at the early phase (assuming no depletion in susceptible fraction).</li>
  <li>Evaluating the effectiveness of case detection, contact tracing and quarantine at different strength level.</li>
  <li>Estimating unknown epidemiological parameters (reproductive numbers without the aforemenioned control measures, proportion of pre-symptomatic transmission, and the generation interval) with cluster and serial interval data.</li>
</ol>

An application on quantifying the effectiveness of COVID-19 control measures using this package can be found in [Here](https://jamanetwork.com/journals/jamainternalmedicine/fullarticle/2778395).
 


### Methods

We adopted the structure of the branching process model developed by [Hellewell, et al.](https://www.thelancet.com/journals/langlo/article/PIIS2214-109X(20)30074-7/fulltext), which in essence consists
of two components: 1) the branching process, and 2) the intervention model. The branching process simulates the
disease transmission dynamics in the early phase of an outbreak. The intervention model considers the effects of
various case-based interventions for COVID-19 control. Specifically, case detection, contact tracing, and quarantine of close contacts are considered.

For the parameter estimation, we applied a Sequential Monte Carlo (SMC) algorithm similar to [Peak, et al.](https://www.pnas.org/content/114/15/4023) since the likelihood of this model is intractable. 
Since different data informs different parts of the model, model fitting will proceed in two stages if both observed serial intervals and cluster size distribution of stuttered transmission chains were presented.
The reproduction number under no interventions (detection, contact tracing and quarantine) and the dispersion parameter  were
estimated by fitting the model to observed cluster size distribution. And the  probability of pre-symptomatic transmission and the standard deviation of the generation interval 
were estimated by fitting the model to the observed serial interval distribution. 

Detailed methods and the parameter for the sampling algorithm can be found [Here.](https://cdn.jamanetwork.com/ama/content_public/journal/intemed/938723/ioi210018supp1_prod_1625245891.34446.pdf?Expires=1641323294&Signature=HsQtyRd33mOSYfaUulXutglY6NxeTlsuI4GCChnNhcxtcxTEfMy1ddALqQxF4G7g1iYZ~f6LjdjVJGUDuEHhgVh2zP4LKfHeIOIkuVKsJqVcuKuNHlGRsj1hWEK4MkkeYx8hFnVaWsVCEcNJiikyfXGirvKii0X9SBV0JlQT7BcREyQC7m-y0-AUA8n09van-8pfFtqhTfA530uSEsvgB2bGMdRDkF8KJPD20VeczmB0XTUcDBbwMA4evYBF5jZGeneAvMnyEpekx5TwwrQBlnaN4hkmK-4Us32yfbGHqA-PqD6N-C~zWTa26xGYeBqsgkdSOKM8MUCRTJzXMOtt4g__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA)



