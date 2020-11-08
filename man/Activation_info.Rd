\name{Activation_info}
\alias{Activation_info}
\docType{data}
\title{ Activation function list
%%   ~~ data name/kind ... ~~
}
\description{This data list the example of the activation function in RNN-based model
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("Activation_info")}
\format{
  A data frame with 277 observations on the following 6 variables.
  \describe{
    \item{\code{X}}{a numeric vector}
    \item{\code{Number}}{a numeric vector}
    \item{\code{Activation_name}}{a factor with levels \code{Act_ART1_NC} \code{Act_ART2_Identity} \code{Act_ART2_NormIP} \code{Act_ART2_NormP} \code{Act_ART2_NormV} \code{Act_ART2_NormW} \code{Act_ART2_Rec} \code{Act_ART2_Rst} \code{Act_ARTMAP_DRho} \code{Act_ARTMAP_NCa} \code{Act_ARTMAP_NCb} \code{Act_at_least_1} \code{Act_at_least_2} \code{Act_at_most_0} \code{Act_BAM} \code{Act_BSB} \code{Act_CC_Thresh} \code{Act_Component} \code{Act_Elliott} \code{Act_Euclid} \code{Act_exactly_1} \code{Act_Exponential} \code{Act_HystStep} \code{Act_Identity} \code{Act_IdentityPlusBias} \code{Act_less_than_0} \code{Act_Logistic} \code{Act_LogisticTbl} \code{Act_LogSym} \code{Act_MinOutPlusWeight} \code{Act_Perceptron} \code{Act_Product} \code{Act_RBF_Gaussian} \code{Act_RBF_MultiQuadratic} \code{Act_RBF_ThinPlateSpline} \code{Act_RM} \code{Act_Signum} \code{Act_Signum0} \code{Act_Sinus} \code{Act_Softmax} \code{Act_StepFunc} \code{Act_TACOMA} \code{Act_TanH} \code{Act_TanHPlusBias} \code{Act_TanH_Xdiv2} \code{Act_TD_Elliott} \code{Act_TD_Logistic} \code{ART1} \code{ART1_Stable} \code{ART1_Synchronous} \code{ART1_Weights} \code{ART2} \code{ART2_Stable} \code{ART2_Synchronous} \code{ART2_Weights} \code{ARTMAP} \code{ARTMAP_Stable} \code{ARTMAP_Synchronous} \code{ARTMAP_Weights} \code{Auto_Synchronous} \code{BackPercolation} \code{BackpropBatch} \code{BackpropChunk} \code{BackpropClassJogChunk} \code{BackpropJogChunk} \code{BackpropMomentum} \code{BackpropWeightDecay} \code{BAM_Order} \code{BBPTT} \code{Binary} \code{BPTT} \code{BPTT_Order} \code{CC} \code{CC_Order} \code{CC_Weights} \code{Clip} \code{ClippHebb} \code{Counterpropagation} \code{CounterPropagation} \code{CPN_Rand_Pat} \code{CPN_Weights_v3.2} \code{CPN_Weights_v3.3} \code{DLVQ_Weights} \code{Dynamic_LVQ} \code{ENZO_noinit} \code{ENZO_prop} \code{Hebb} \code{Hebb_Fixed_Act} \code{Hebbian} \code{Hopfield_Fixed_Act} \code{Hopfield_Synchronous} \code{Inverse} \code{JE_BP} \code{JE_BP_Momentum} \code{JE_Order} \code{JE_Quickprop} \code{JE_Rprop} \code{JE_Special} \code{JE_Weights} \code{Kohonen} \code{Kohonen_Const} \code{Kohonen_Order} \code{Kohonen_Rand_Pat} \code{Kohonen_Weights_v3.2} \code{LinearScale} \code{Logistic_notInhibit} \code{MagPruning} \code{Monte-Carlo} \code{Noncontributing_Units} \code{None} \code{Norm} \code{OptimalBrainDamage} \code{OptimalBrainSurgeon} \code{Out_ART2_Noise_ContDiff} \code{Out_ART2_Noise_PLin} \code{Out_Clip_01} \code{Out_Clip_11} \code{Out_Identity} \code{Out_Threshold05} \code{PruningFeedForward} \code{PseudoInv} \code{QPTT} \code{Quickprop} \code{RadialBasisLearning} \code{Randomize_Weights} \code{Random_Order} \code{Random_Permutation} \code{Random_Weights_Perc} \code{RBF-DDA} \code{RBF_Weights} \code{RBF_Weights_Kohonen} \code{RBF_Weights_Redo} \code{RM_delta} \code{RM_Random_Weights} \code{Rprop} \code{RpropMAP} \code{SCG} \code{Serial_Order} \code{Sim_Ann_SS} \code{Sim_Ann_WTA} \code{Sim_Ann_WWTA} \code{Site_at_least_1} \code{Site_at_least_2} \code{Site_at_most_0} \code{Site_Max} \code{Site_Min} \code{Site_Pi} \code{Site_Produkt} \code{Site_Reciprocal} \code{Site_WeightedSum} \code{Skeletonization} \code{Std_Backpropagation} \code{Synchonous_Order} \code{TACOMA} \code{Threshold} \code{TimeDelayBackprop} \code{TimeDelay_Order} \code{Topological_Order}}
    \item{\code{type}}{a numeric vector}
    \item{\code{Number_of_inParams}}{a numeric vector}
    \item{\code{Number_of_outParams}}{a numeric vector}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(Activation_info)
## maybe str(Activation_info) ; plot(Activation_info) ...
}
\keyword{datasets}
