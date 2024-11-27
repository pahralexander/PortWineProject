using JSON
using Distributions
using Dates
using Profile
using Base.Threads
using Random
#using PyPlot
using StatsBase
using SparseArrays
using Combinatorics
#----------------------------------#
function integratedBlendingAndPurchasing()
	for nP in [2]
		for sF in [25/6.0]
			for sprRev in [3]
				for pU in [true]
					for oDP in [0.25]
						priceRange = pU ? [20] : ["na"]
						for pS in priceRange
						resultsDict = Dict()
						pUString = pU ? "PU" : "NoPU"
						decayProbs = Dict{Int32,SparseVector{Float64,Int32}}(); dSupermarket = Dict{Int32,Float64}(); #dStates = Dict{Int32,Array{Int32,1}}();
						currentInstance = string(string(@__DIR__),"/Results/",nP,"nP_",sF,"sF_",pUString,"pU_",sprRev,"sprRev_",oDP,"oDP_",pS,"pS_MDP_approximations.json")
						for approx in ["MDP"]

								# nP = 2
								# sF = 3
								# sprRev = 3
								# oDP = 0.2
								# pS = 30

								start = now()

								# open(currentInstance,"r") do f
								# 	resultsDict = JSON.parse(f)
								# end

								ages = [i for i in 1:6]
								numAges = ages[end]

								scalingFactor = sF

								agesLarge = [i for i in 1:25]
								numAgesLarge = agesLarge[end]

								largeSpread = 10

								ageBuckets = [1,1,1,1,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6]

								#approximation case
								BAp = approx == "BAp" ? true : false
								SAp = approx == "SAp" ? true : false
								PrAp = approx == "PrAp" ? true : false
								PuAp = approx == "PuAp" ? true : false
								PuSAp = approx == "PuSAp" ? true : false
								PuPrAp = approx == "PuPrAp" ? true : false
								TAp = approx == "TAp" ? true : false


								#initialize data from original large scale case
								supermarketStart = 30/scalingFactor
								supermarketStep = 30/scalingFactor
								maxAgeSupermarket = 3*scalingFactor
								holdingCostsLarge = 30/scalingFactor
								supermarketMax = maxAgeSupermarket*supermarketStep
								supermarketCont = [min(supermarketStart+(i-1)*supermarketStep,supermarketMax) for i in 1:numAgesLarge]
								stockValue = [holdingCostsLarge * i for i in 1:numAgesLarge]
								supermarketContributionLarge = supermarketCont-stockValue

								#define decay probability in large scale case
								qWeib = 0.875
								βWeib = 0.8
								overallDecay = oDP
								decayProbabilityLarge = zeros(numAgesLarge)
								approxdecayProb = zeros(numAgesLarge)
								for k in 1:agesLarge[end]
									approxdecayProb[k] = qWeib^((k-1)^(βWeib))-qWeib^((k)^(βWeib))
								end
								for k in 1:agesLarge[end]
								   	decayProbabilityLarge[k] = (approxdecayProb[k]/sum(values(approxdecayProb)))*overallDecay
								end
								println(decayProbabilityLarge, sum(decayProbabilityLarge))

								#map supermarket contribution and decay proabability from large scale case to base case
								currAgeClass = 1
								decayProbability = zeros(numAges)
								supermarketContribution = zeros(numAges)
								currDecayVector = []
								currSMVector = []
								for i in 1:numAgesLarge
									nextAgeClass = ageBuckets[i]
									if nextAgeClass == currAgeClass
										append!(currDecayVector,1-decayProbabilityLarge[i])
										append!(currSMVector,supermarketContributionLarge[i])
									else
										decayProbability[currAgeClass] = 1-prod(currDecayVector)
										supermarketContribution[currAgeClass] = sum(currSMVector[k]*(currDecayVector[k]/sum(currDecayVector)) for k in 1:length(currDecayVector))
										currAgeClass = nextAgeClass
										currDecayVector = [1-decayProbabilityLarge[i]]
										currSMVector = [supermarketContributionLarge[i]]
									end
									if i == numAgesLarge
										decayProbability[currAgeClass] = 1-prod(currDecayVector)
										supermarketContribution[currAgeClass] = sum(currSMVector[k]*(currDecayVector[k]/sum(currDecayVector)) for k in 1:length(currDecayVector))
									end
								end
								println(supermarketContribution)
								println(decayProbability, sum(decayProbability))


								#initialize Wine classes
								nProducts = nP
								wineClasses = [i for i in 1:nProducts]
								targetAges = nProducts > 2 ? [3,4,5] : [3,5]
								detDemand = [4,4]
								totalDemand = sum(detDemand)

								#initialize Inventory levels
								nInventoryLevels = 13
								nPriceLevels = pU ? 9 : 1
								stateMultipliers = [nInventoryLevels^(numAges-i) for i in 1:numAges]

								#data representation for saving instance data
								contrExp = sqrt(sprRev)
								contrStart = sprRev == 4 ? 800/3 : 1000/3
								brandContribution = [contrStart*contrExp^(targetAges[i]-targetAges[1]) for i in 1:nProducts]
								holdingCosts = holdingCostsLarge * scalingFactor

								HtoM = round(brandContribution[end]/brandContribution[1], digits=1)
								HtoS = round(brandContribution[end]/supermarketContribution[end], digits=1)
								HtoHC = round(brandContribution[end]/holdingCosts, digits=1)

								muH = (nPriceLevels-1)/2
								sigmaH = 2.04058

								#wine costs (=purchasing costs) of each respective age class
								wineCostMean = 150
								costStep = pS
								if pU
									harvestCosts = nPriceLevels%2 == 0 ? [] : [wineCostMean]
									for i in 1:floor(nPriceLevels/2)
										append!(harvestCosts,wineCostMean - (i-(((nPriceLevels+1)%2)/2)) * costStep)
										append!(harvestCosts,wineCostMean + (i-(((nPriceLevels+1)%2)/2)) * costStep)
									end
								else
									harvestCosts = [wineCostMean]
								end
								sort!(harvestCosts)
								HtoPC = round(brandContribution[end]/mean(harvestCosts), digits=1)
								println(harvestCosts)

								#maximally allowed spread in blending
								maxSpread = 2

								#indicator for number of harvest price scenarios
								minYield = 0


								if PuAp || TAp || BAp || PrAp || SAp || PuSAp ||PuPrAp
									simulationInstance = ""
									rootComponents=splitpath(@__DIR__)
									for i in 1:length(rootComponents)
										if i==1
											simulationInstance *= rootComponents[i][1:2]*"/"
										elseif i==length(rootComponents)
											simulationInstance *= rootComponents[i]
										else
											simulationInstance *= rootComponents[i]*"/"
										end
									end
									pResults = Dict()
									if TAp
										MDPinstance = simulationInstance * string("/Results/",nP,"nP_",sF,"sF_",pUString,"pU_",sprRev,"sprRev_",oDP,"oDP_",pS,"pS_MDP.json")
										open(MDPinstance, "r") do f
											optResults = JSON.parse(f)
										end
										optimalAvgReward = optResults["averageReward"]
										resultsDict["MDP"] = optimalAvgReward; optResults = nothing
									end
									simulationInstance = simulationInstance * string("/Analysis/Simulation_Aggregated/",nP,"nP_",sF,"sF_",pUString,"pU_",sprRev,"sprRev_",oDP,"oDP_",pS,"pS_Simulation.json")
									open(simulationInstance, "r") do f
										pResults = JSON.parse(f)
									end
									purchasingResults = pResults["priceStockPurchasing"]
									purchasing3Results = pResults["priceStock3Purchasing"]
									rule2D = Dict(p=>pResults["2DRule"][string(p)] for p in 0:nPriceLevels-1)

									for p in keys(purchasingResults) for k in keys(purchasingResults[p]) if purchasingResults[p][k] == nothing purchasingResults[p][k] = -1 end end end
									for p in keys(purchasing3Results) for k in keys(purchasing3Results[p]) if purchasing3Results[p][k] == nothing purchasing3Results[p][k] = -1 end end end
									optimalBlendingInv = pResults["ageStructureBlending"]
									optimalPurchasingInv = vcat([pResults["averagePurchasing"]],pResults["ageStructurePurchasing"])
									optimalPostBlendingInv = pResults["ageStructurePurchasing"]
								end

								# helper data for defining actions
								# define maximum sales for each age class
								maxSales = Dict(i=>0 for i in ages)
								for i in ages
									if i >= targetAges[1]
										maxSales[i]+=sum(detDemand[j] for j in wineClasses if i>=targetAges[j])
									end
									if i < targetAges[end]
										nextProd = findfirst(j -> j > i,targetAges)
										for j in nextProd:nProducts
											for k in (detDemand[j]-1):-1:1
												if k*i+(detDemand[j]-k)*min(maxSpread+i,numAges) >= targetAges[j]*detDemand[j]
													maxSales[i] += k
													break
												end
											end
										end
									end
								end
								#println("maxSales ",maxSales)

								#define minimum sum of all ages for different demand scenarios

								demandPatterns = [[i] for i in 0:detDemand[1]]
								for p in 2:nProducts
									currLength = length(demandPatterns)
									for i in 0:detDemand[p]
										for k in 1:currLength
											if i == detDemand[p] && k == currLength
												toExtend = copy(demandPatterns[k])
												append!(toExtend,i)
												deleteat!(demandPatterns,1:k)
												demandPatterns = vcat(demandPatterns,[toExtend])
											else
												toExtend = copy(demandPatterns[k])
												append!(toExtend,i)
												demandPatterns = vcat(demandPatterns,[toExtend])
											end
										end
									end
								end
								println(demandPatterns)
								#define possible age patterns for different sales quantities
								possiblePatterns = Dict(i=>[] for i in 0:totalDemand)
								for d in demandPatterns
									append!(possiblePatterns[sum(d)],[d])
								end
								minimumAgeSum=Dict(i=>Inf for i in 0:totalDemand)
								for d in demandPatterns
									if d' * targetAges < minimumAgeSum[sum(d)]
										minimumAgeSum[sum(d)] = d' * targetAges
									end
								end
								println("possiblePatterns ",possiblePatterns)

								#define decay probabilities (use discrete Weibull distribution as in Nakagawa and Osaki (1975)) + Bernoulli-process resulting partial decay probabilities
								# decayProbability = zeros(numAges)
								# approxdecayProb = zeros(numAges)
								# for k in 1:ages[end]
								# 	approxdecayProb[k] = qWeib^((k-1)^(βWeib))-qWeib^((k)^(βWeib))
								# end
								# for k in 1:ages[end]
								#    	decayProbability[k] = round((approxdecayProb[k]/sum(values(approxdecayProb)))*overallDecay,digits=3)
								# end

								println("DECAY PROBABAILITY:\n",decayProbability, sum(values(decayProbability)))
								combinedDecay = Dict(s=>zeros(numAges,s+1) for s in 1:nInventoryLevels-1)
								for s in 1:nInventoryLevels-1
									for d in 0:s
										for a in 1:numAges
											combinedDecay[s][a,d+1] = binomial(s,d)*(decayProbability[a]^d)*((1-decayProbability[a])^(s-d))
										end
									end
								end
								final = now()-start
								println("Time to set up data and probabilities: ",final)

								#define harvest yield probabilities for price definition (truncating a discrete normal distribution between 0 and 6)
								if pU
									lambda = exp((2*muH-1)/(2*(sigmaH^2))); qNorm = exp(-1/(sigmaH^2))
									pmfDivisor=sum((lambda^j)*(qNorm^(j*(j-1)/2)) for j in -200:200)
									untruncYield = Dict(i=>((lambda^i)*(qNorm^(i*(i-1)/2)))/pmfDivisor for i in 0:nPriceLevels-1)
									yieldProbability = Dict(i=>untruncYield[i]/sum(values(untruncYield)) for i in 0:nPriceLevels-1)
									println(yieldProbability, sum(values(yieldProbability)))
								else
									yieldProbability = Dict(0=>1)
								end

								#MODEL
								#state space, action space, and decays
								#define states

								start = now()
								allStates = defineStateSpace(0,numAges,nInventoryLevels)
								stateSpaceSize = length(allStates[:,1])
								stateSpace = Dict{Int32,Array{Int8,1}}(); sizehint!(stateSpace,stateSpaceSize)
								stateSpace = Dict(((allStates[s,:])' * stateMultipliers)+1=>allStates[s,:] for s in 1:stateSpaceSize)
								purchasingStateSpace = Dict()
								reversedStateSpace = Dict(v=>k for (k,v) in stateSpace)
								println("Number of states: ", stateSpaceSize)
								final=now()-start
								println("Time to create state space: ", final)

								# println(isReasonable([0,0,0,0,1],numAges,minimumAgeSum,possiblePatterns,targetAges,detDemand,maxSpread,wineClasses))
								# error("stop here")

								#define blending actions
								start=now()
								allActions = defineActionSpace(1,maxSales,numAges,totalDemand,minimumAgeSum,possiblePatterns,targetAges,detDemand,maxSpread,wineClasses,stateSpace)
								actionSpaceSize = length(allActions[:,1])
								println("Number of actions: ", actionSpaceSize)
								blendingActionSpace = Dict(a=>allActions[a,1:numAges] for a in 1:actionSpaceSize)
								stateEncodedBlendings = Dict(a=>(blendingActionSpace[a])' * stateMultipliers +1 for a in 1:actionSpaceSize)
								demandFulfillment = Dict(a=>allActions[a,numAges+1:numAges+length(wineClasses)] for a in 1:actionSpaceSize)
								prodSpecActions = Dict(a=>Dict(w=>allActions[a,numAges*w+length(wineClasses)+1:numAges*(w+1)+length(wineClasses)] for w in wineClasses) for a in 1:actionSpaceSize)

								# for a in 1:actionSpaceSize
								# 	println(a, ": ", blendingActionSpace[a])
								# end

								averageAge = Dict(a => Dict(w => ((prodSpecActions[a][w])' * ages)/max(demandFulfillment[a][w],1) for w in wineClasses) for a in 1:actionSpaceSize)
								blendingProduct = Dict(a => Dict(w => getBlending(prodSpecActions[a][w]) for w in wineClasses) for a in 1:actionSpaceSize)
								fulfillProduct = Dict(a => Dict(w => demandFulfillment[a][w] for w in wineClasses) for a in 1:actionSpaceSize)

								actionContribution = Dict(a => sum(fulfillProduct[a][w]*(brandContribution[w] - ((averageAge[a][w]-targetAges[w])*holdingCosts)) for w in wineClasses) for a in 1:actionSpaceSize); allActions=nothing

								#define purchasing actions
								purchasingActionSpace = [i for i in 0:nInventoryLevels-1]
								counter = 0
								for p in 0:nPriceLevels-1
									for i in 1:stateMultipliers[1]
										counter += 1
										purchasingStateSpace[counter] = vcat([p],stateSpace[i][2:end])
									end
								end
								purchasingStateSpaceSize = counter
								dSupermarketPurchase = zeros(length(purchasingActionSpace))
								age1Probs = zeros(length(purchasingActionSpace),length(purchasingActionSpace))
								for p in purchasingActionSpace
									if p == 0
										dSupermarketPurchase[p+1] = 0
										age1Probs[1,1] = 1.0
									else
										dSupermarketPurchase[p+1] = foldl(+,collect(combinedDecay[p][1,d+1]*(d*supermarketContribution[1]) for d in 0:p))
										for r in 0:p
											age1Probs[p+1,r+1] = combinedDecay[p][1,p-r+1]
										end
									end
								end
								println(purchasingStateSpaceSize)
								statePrice = [harvestCosts[purchasingStateSpace[s][1]+1] for s in 1:purchasingStateSpaceSize]
								
								#define decays based on pre-purchase states
								if approx in ["PuAp","BAp", "MDP"]
									start = now()
									transStates = [allStates[i,2:end] for i in 1:stateMultipliers[1]]
									#sizehint!(decayProbs,stateMultipliers[1]); #sizehint!(dSupermarket,stateMultipliers[1]); #sizehint!(dStates,stateSpaceSize);
									println("Decay transition generation started!")
									decayProbs = Dict(i => spzeros(stateMultipliers[1]) for i in 1:stateMultipliers[1])
									#decayContr = spzeros(stateMultipliers[1])
									for t in 1:stateMultipliers[1]
										dSupermarket[t] = 0.0
										input = copy(transStates[t])
										for tTrans in 1:t
											transState = transStates[tTrans]
											if isTransState(input, transState)
												decayProb = 1.0
												decayContr = 0.0
												for i in 1:numAges-1
													if input[i] > 0
														decay = input[i] - transState[i]
														decayProb *= combinedDecay[input[i]][i+1, decay + 1]
														decayContr += decay * supermarketContribution[i+1] 
													end
												end
												dSupermarket[t] += decayContr * decayProb
												decayProbs[t][tTrans] = decayProb
											end
										end
										if t % 500 == 0
											println(t, " states done in: ", now() - start)
										end
										#if !isapprox(sum(decayProbs[t]),1,atol=0.001) error("probs don't add up to 1, but ", sum(decayProbs[t]), " at a length of ", length(decayProbs[t]), " at t = ", t); end
									end
								end

								final = now()-start; tempTrans = nothing; transInputs = nothing; reversedStateSpace = nothing; tStates = nothing; nTrans=nothing
								println("Time to create decay transition space: ", final)

								#define action-state rewards for initial model + applicable actions for each state
								stateBlendingActions = Dict{Int32,Array{Int16,1}}(); statePurchasingInfo = Dict{Int32,Dict{Int16,Float64}}(); stateBlendingInfo=Dict{Int32,Dict{Int16,Tuple{Float64,Int32}}}();
								sizehint!(stateBlendingActions,stateSpaceSize); sizehint!(statePurchasingInfo,purchasingStateSpaceSize); sizehint!(stateBlendingInfo,stateSpaceSize);
								for s in Int32(1):Int32(stateSpaceSize)
									stateBlendingActions[s] = Array{Int16,1}(); stateBlendingInfo[s] = Dict{Int16,Tuple{Float64,Int32}}()
									bestApproxPattern = zeros(nProducts)
									for a in Int16(1):Int16(actionSpaceSize)
										if isCompatible(stateSpace[s],blendingActionSpace[a],totalDemand)
											#actions may only be executed if demand is satisfied to the greatest extent possible
											if BAp || SAp || TAp || PuSAp
												if reverse(demandFulfillment[a]) > bestApproxPattern
													stateBlendingActions[s] = []
													stateBlendingInfo[s] = Dict()
													bestApproxPattern = demandFulfillment[a]
												end
											end
											#get deterministic outcomes from state-action-pair
											newState = stateSpace[s] - blendingActionSpace[a]
											carryOnState = newState[1:end-1]
											#fetch transition information based on interim inventory state
											stateBlendingInfo[s][a] = actionContribution[a], (carryOnState)' * stateMultipliers[2:end] +1 #+ newState[end]*supermarketContribution[end] no outdating contributions
											append!(stateBlendingActions[s],a)
										end
									end
									if isempty(stateBlendingInfo[s])
										println(stateSpace[s])
										error("no applicable action")
									end
									if BAp || PrAp || TAp || PuPrAp
										if BAp || TAp
											actionAlternatives = copy(stateBlendingActions[s])
											#maximize target age adherence from older to younger wine
											for w in wineClasses[end]:-1:1
												salesProduct = fulfillProduct[actionAlternatives[end]][w]
												remAlternatives = []
												maxTargetAdh = numAges - targetAges[w]
												for a in actionAlternatives
													targetAdh = (sum(prodSpecActions[a][w][i]*i for i in ages)/max(1,salesProduct)) - targetAges[w]
													if targetAdh < maxTargetAdh
														maxTargetAdh = targetAdh
														remAlternatives = [a]
													elseif targetAdh == maxTargetAdh
														append!(remAlternatives,a)
													end
												end
												actionAlternatives = copy(remAlternatives)
											end
											for w in wineClasses[end]:-1:1
												minSpread = numAges
												remAlternatives = []
												if demandFulfillment[actionAlternatives[1]][w] > 0
													for a in actionAlternatives
														spread = findlast(x -> x>0, prodSpecActions[a][w]) - findfirst(x -> x>0, prodSpecActions[a][w])
														if spread < minSpread
															minSpread = spread
															remAlternatives = [a]
														elseif spread == minSpread
															append!(remAlternatives,a)
														end
													end
													actionAlternatives = copy(remAlternatives)
												end
											end
											minDifference = Inf
											finalAlternatives = copy(actionAlternatives)
											# optimalPurchasingInv = []
											# for a in actionAlternatives
											# 	nextState = stateSpace[s]-blendingActionSpace[a]
											# 	differenceFromOptimal = sum(map(x -> abs(x)^2,nextState - optimalPostBlendingInv))
											# 	if differenceFromOptimal < minDifference
											# 		minDifference = differenceFromOptimal
											# 		finalAlternatives = [a]
											# 	end
											# end
											if length(finalAlternatives)>1
												error("too many alternatives in state ", string(stateSpace[s]), string([blendingActionSpace[a] for a in remAlternatives4]))
											elseif finalAlternatives == []
												error("no alternatives in state ", string(stateSpace[s]))
											end
											stateBlendingActions[s] = finalAlternatives

											stateBlendingInfo[s] = Dict()
											stateBlendingInfo[s][finalAlternatives[1]] = actionContribution[finalAlternatives[1]] + (stateSpace[s]-blendingActionSpace[finalAlternatives[1]])[end]*supermarketContribution[end], ((stateSpace[s]-blendingActionSpace[finalAlternatives[1]])[1:end-1])' * stateMultipliers[2:end] +1
										else
											dAlternatives = Dict(d=>[] for d in demandPatterns)
											stateBlendingInfo[s] = Dict()
											for a in stateBlendingActions[s]
												append!(dAlternatives[demandFulfillment[a]],a)
											end
											stateBlendingActions[s] = []
											for d in demandPatterns
												if dAlternatives[d] != []
													actionAlternatives = copy(dAlternatives[d])
													#maximize target age adherence from older to younger wine
													for w in wineClasses[end]:-1:1
														maxTargetAdh = nInventoryLevels
														remAlternatives = []
														for a in actionAlternatives
															targetAdh = (sum(prodSpecActions[a][w][i]*i for i in ages)/max(1,sum(d))) - targetAges[w]
															if targetAdh < maxTargetAdh
																maxTargetAdh = targetAdh
																remAlternatives = [a]
															elseif targetAdh == maxTargetAdh
																append!(remAlternatives,a)
															end
														end
														actionAlternatives = copy(remAlternatives)
													end
													for w in wineClasses[end]:-1:1
														minSpread = numAges
														remAlternatives = []
														if demandFulfillment[actionAlternatives[1]][w] > 0
															for a in actionAlternatives
																spread = findlast(x -> x>0, prodSpecActions[a][w]) - findfirst(x -> x>0, prodSpecActions[a][w])
																if spread < minSpread
																	minSpread = spread
																	remAlternatives = [a]
																elseif spread == minSpread
																	append!(remAlternatives,a)
																end
															end
															actionAlternatives = copy(remAlternatives)
														end
													end
													minDifference = Inf
													finalAlternatives = copy(actionAlternatives)
													if length(finalAlternatives) > 1
														error("too many alternatives")
													end
													# optimalPurchasingInv = []
													# for a in actionAlternatives
													# 	nextState = stateSpace[s]-blendingActionSpace[a]
													# 	differenceFromOptimal = sum(map(x -> abs(x)^2,nextState - optimalPostBlendingInv))
													# 	if differenceFromOptimal < minDifference
													# 		minDifference = differenceFromOptimal
													# 		finalAlternatives = [a]
													# 	end
													# end
													append!(stateBlendingActions[s],finalAlternatives[1])
													stateBlendingInfo[s][finalAlternatives[1]] = actionContribution[finalAlternatives[1]] + (stateSpace[s]-blendingActionSpace[finalAlternatives[1]])[end]*supermarketContribution[end], ((stateSpace[s]-blendingActionSpace[finalAlternatives[1]])[1:end-1])' * stateMultipliers[2:end] +1
												end
											end
										end
									end
								end
								for s in 1:purchasingStateSpaceSize
									statePurchasingInfo[s] = Dict{Int16,Float64}();
									if PuAp || TAp || PuSAp || PuPrAp
										currPrice = purchasingStateSpace[s][1]
										invT = sum(purchasingStateSpace[s][2:end])
										inv3 = sum(purchasingStateSpace[s][2:3])
										# pRuleT = Int(round(purchasingResults[string(purchasingStateSpace[s][1])][string(sum(purchasingStateSpace[s][2:end]))], digits=0))
										# pRule3 = Int(round(purchasing3Results[string(purchasingStateSpace[s][1])][string(sum(purchasingStateSpace[s][2:3]))], digits=0))
										# if pRuleT < 0
										# 	if sum(purchasingStateSpace[s][2:end]) > sum(optimalPostBlendingInv)
										# 		pRule = 0
										# 	else
										# 		pRule = nInventoryLevels-1
										# 	end
										# else
										# 	pRule = pRuleT
										# end

										pRule = min(max(0,rule2D[currPrice][1] - invT, rule2D[currPrice][2] - inv3),nInventoryLevels-1)

										prePurchaseState = (purchasingStateSpace[s][2:end])' * stateMultipliers[2:end] + 1
										statePurchasingInfo[s][pRule] = -statePrice[s]*pRule - holdingCosts*(sum(purchasingStateSpace[s][2:end])+pRule) + dSupermarketPurchase[pRule+1] + dSupermarket[prePurchaseState]
									else
										for p in purchasingActionSpace
											prePurchaseState = (purchasingStateSpace[s][2:end])' * stateMultipliers[2:end] + 1
											statePurchasingInfo[s][p] = -statePrice[s]*p - holdingCosts*(sum(purchasingStateSpace[s][2:end])+p) + dSupermarketPurchase[p+1] + dSupermarket[prePurchaseState]
										end
									end
								end
								tSupermarketInit = nothing;
								countActions = 0
								for s in 1:stateSpaceSize
									countActions += length(stateBlendingActions[s])
								end
								println("Avg number of actions per state: ", countActions/stateSpaceSize)
								final = now() - start; demandFulfillment=nothing
								println("Time to create state specific actions space: ", final)

								if TAp
									currentInstance = string(string(@__DIR__),"/Results/",nP,"nP_",sF,"sF_",pUString,"pU_",sprRev,"sprRev_",oDP,"oDP_",pS,"pS_TAp.json")
									TApPolicyB = Dict(s => stateBlendingActions[s][1] for s in 1:stateSpaceSize)
									TApPolicyP = Dict(s => collect(keys(statePurchasingInfo[s]))[1] for s in 1:purchasingStateSpaceSize)
									TAPPolicies = Dict("blending" => TApPolicyB, "purchasing"=>TApPolicyP)
									open(currentInstance, "w") do f
										JSON.print(f,TAPPolicies,4)
									end
								end

								#vALUE ITERATION ALGORITHM
								optimalBlending, optimalPurchasing, averageReward = valueIterationIntegrated(stateBlendingInfo,statePurchasingInfo,decayProbs,stateSpaceSize,purchasingStateSpaceSize,yieldProbability,stateMultipliers,age1Probs)
								#
								#PRINT OUT RESULTS
								if SAp || PrAp || BAp || PuAp || TAp || PuSAp || PuPrAp
									currentInstance = string(string(@__DIR__),"/Results/",nP,"nP_",sF,"sF_",pUString,"pU_",sprRev,"sprRev_",oDP,"oDP_",pS,"pS_MDP_approximations.json")
									resultsDict[approx] = averageReward
								else
									currentInstance = string(string(@__DIR__),"/Results/",nP,"nP_",sF,"sF_",pUString,"pU_",sprRev,"sprRev_",oDP,"oDP_",pS,"pS_MDP.json")
									resultsDict = Dict("stateSpace"=>stateSpace,"purchasingStateSpace"=>purchasingStateSpace,"blendingActionSpace"=>blendingActionSpace,"fulfillProduct"=>fulfillProduct,"optimalBlending"=>optimalBlending,"optimalPurchasing"=>optimalPurchasing,"averageReward"=>averageReward,"stateBlendingActions"=>stateBlendingActions,"blendingProduct"=>blendingProduct,"prodSpecActions"=>prodSpecActions)
									open(currentInstance,"w") do f
										JSON.print(f,resultsDict,4)
									end
								end
							end
							# if haskey(resultsDict,"TAp")
							# 	open(currentInstance,"w") do f
							# 		JSON.print(f,resultsDict,4)
							# 	end
							# end
						end
					end
				end
			end
		end
	end
end
#----------------------------------------#
# value iteration algorithm for integrated blending and purchasing model
function valueIterationIntegrated(stateBlendingInfo,statePurchasingInfo,decayProbs,stateSpaceSize,purchasingStateSpaceSize,yieldProbability,stateMultipliers,age1Probs)
	start=now()
	convergenceBlending = Inf; convergencePurchasing = Inf;
	epsilon = 0.001
	ValueBlending = zeros(Float64,2,stateSpaceSize)
	ValuePurchasing = zeros(Float64,2,purchasingStateSpaceSize)
	BlendingPolicy = zeros(Int32,stateSpaceSize)
	PurchasingPolicy = zeros(Int32,purchasingStateSpaceSize)
	iteration = 0
	numCarryOnStates = stateMultipliers[1]
	prePurchaseState = [((s-1)%stateMultipliers[1])+1 for s in 1:stateSpaceSize]
	statePurchasing = [Int(floor((s-1)/stateMultipliers[1])) for s in 1:stateSpaceSize]
	purchasingActions = [i for i in 0:length(age1Probs[:,1])-1]
	purchasingSections = [i*stateMultipliers[1]+1:i*stateMultipliers[1]+stateMultipliers[1] for i in 0:length(purchasingActions)-1]
	averageReward = -Inf

	while convergenceBlending > epsilon || convergencePurchasing > epsilon
		iteration += 1
		println("Iteration: ", iteration)

		previousPurchasingValues = view(ValuePurchasing,(iteration-1)%2+1,:)
		expectedPreviousPurchasingValues = Dict(); sizehint!(expectedPreviousPurchasingValues,numCarryOnStates)
		for t in 1:numCarryOnStates
			expectedPreviousPurchasingValues[t] = sum(yieldProbability[y] * previousPurchasingValues[t+y*stateMultipliers[1]] for y in 0:length(yieldProbability)-1)
		end
		#compute actions for blending states
		@inbounds for s in Int32(1):Int32(stateSpaceSize)
			ValueBlending[iteration%2+1,s] = -Inf
			for (a,i) in stateBlendingInfo[s]
				ActionValue = i[1] + expectedPreviousPurchasingValues[i[2]]
				if ActionValue > ValueBlending[iteration%2+1,s]
					ValueBlending[iteration%2+1,s] = ActionValue
					BlendingPolicy[s] = a
				end
			end
		end
		println(minimum(ValueBlending[iteration%2+1,:]))
		println("Blending states done after ", now()-start)
		previousBlendingValues = [foldl(+,[view(ValueBlending,iteration%2+1,purchasingSections[r+1])*age1Probs[p+1,r+1] for r in 0:p]) for p in purchasingActions]
		println("fold left for purchasing actions done after ", now() - start)
		expectedPreviousBlendingValues = [(previousBlendingValues[statePurchasing[s]+1])' * decayProbs[prePurchaseState[s]] for s in 1:stateSpaceSize]

		#compute actions for purchasing states
		@inbounds for s in Int32(1):Int32(purchasingStateSpaceSize)
			ValuePurchasing[iteration%2+1,s] = -Inf
			for (p,v) in statePurchasingInfo[s]
				ActionValue = v + expectedPreviousBlendingValues[p*stateMultipliers[1] + prePurchaseState[s]]
				if ActionValue > ValuePurchasing[iteration%2+1,s]
					ValuePurchasing[iteration%2+1,s] = ActionValue
					PurchasingPolicy[s] = p
				end
			end
		end
		println(minimum(ValuePurchasing[iteration%2+1,:]))
		println("Purchasing states done after ", now()-start)


		MBlending = maximum(view(ValueBlending,iteration%2+1,:) - view(ValueBlending,(iteration-1)%2+1,:))
		mBlending = minimum(view(ValueBlending,iteration%2+1,:) - view(ValueBlending,(iteration-1)%2+1,:))
		convergenceBlending = (MBlending-mBlending)/mBlending

		MPurchasing = maximum(view(ValuePurchasing,iteration%2+1,:) - view(ValuePurchasing,(iteration-1)%2+1,:))
		mPurchasing = minimum(view(ValuePurchasing,iteration%2+1,:) - view(ValuePurchasing,(iteration-1)%2+1,:))
		convergencePurchasing = abs((MPurchasing - mPurchasing)/mPurchasing)
	    println("convergenceBlending: ",convergenceBlending)
		println(mBlending)
		println("convergencePurchasing: ",convergencePurchasing)
		println(mPurchasing)
		println("Time: ",now()-start)
		averageReward = mPurchasing
	end
	final = now() - start
	println("Solved in ", final)
	println("Completed in ", iteration, " iterations")

	return Dict(i => BlendingPolicy[i] for i in 1:stateSpaceSize), Dict(i => PurchasingPolicy[i] for i in 1:purchasingStateSpaceSize), averageReward

end
#----------------------------------------#
# recursive function for the generation of the state space
function defineStateSpace(seed, numAges, nInventoryLevels, stateVector = Array{Int8,1}(zeros(numAges)), q=Int8(-1))
	sV = copy(stateVector)
	pos = copy(seed)
	if pos > 0
		sV[pos] = q
	end
	if pos == numAges
		return transpose(sV)
	end
	currentElement = pos+1
	stateSpace = Array{Int8}(undef,0,numAges)
	for q in 0:nInventoryLevels-1
		newVectors = defineStateSpace(pos+1,numAges,nInventoryLevels,sV,Int8(q))
		if newVectors != nothing
			stateSpace = vcat(stateSpace,newVectors)
		end
	end
	return stateSpace
end
#------------------------------------------------#
#recursive function for the configuration of all actions
function defineActionSpace(seed,maxSales,numAges,totalDemand,minimumAgeSum,possiblePatterns,targetAges,detDemand,maxSpread,wineClasses,stateSpace)#,actionVector=Array{Int8,1}(zeros(numAges*(length(wineClasses)+1)+length(wineClasses))))
	pos = copy(seed)
	stateSpaceSize = length(keys(stateSpace))
	actionSpace = Array{Int8}(undef,0,numAges*(length(wineClasses)+1)+length(wineClasses))
	for s in 1:stateSpaceSize
		if sum(stateSpace[s]) <= totalDemand
			actionVector = vcat(copy(stateSpace[s]),zeros((numAges+1)*length(wineClasses)))
			check = isReasonable(copy(stateSpace[s]),numAges,minimumAgeSum,possiblePatterns,targetAges,detDemand,maxSpread,wineClasses)
			if check != nothing
				aV = copy(actionVector)
				for p in 1:length(check[1])
					aV[numAges+p] = check[1][p]
				end
				for i in 1:numAges
					for w in wineClasses
						aV[numAges*w+length(check[1])+i] = check[2][w][i]
					end
				end
				actionSpace = vcat(actionSpace,transpose(aV))
			end
		end
	end
	return actionSpace
end
#------------------------------------------------#
#checks whether a particular action is valid
function isReasonable(actionVector,numAges,minimumAgeSum,possiblePatterns,targetAges,detDemand,maxSpread,wineClasses)
	daV = copy(actionVector)
	salesAmount = sum(daV)
	ageSum=sum(i*daV[i] for i in 1:numAges)
	if ageSum < minimumAgeSum[salesAmount]
		return nothing
	end
	salesAges = getBlending(daV)

	#check if the wine ages are compatible with a possible sales pattern
	for p in reverse(possiblePatterns[salesAmount])
		pOrig = Array{Int64,1}(copy(salesAges))
		prodSpecActions = [zeros(numAges) for w in wineClasses]
		usedAges = [Array{Int64,1}(zeros(p[w])) for w in wineClasses]
		#first fulfill demand for younger age port wines with the minimum average age possible
		for w in wineClasses
			minSum = Inf
			if p[w] > 0
				for k in powerset(pOrig,p[w],p[w])
					if sum(k) >= targetAges[w]*p[w] && sum(k) < minSum
						if p[w] <= 1 || !any(k[i]-k[1] > maxSpread for i in 2:length(k))
							usedAges[w] = k
							minSum = sum(k)
						end
					end
				end
				if usedAges[w] != zeros(p[w])
					for i in usedAges[w]
						remove!(pOrig,i)
					end
				end
			end
		end
		if !any(sum(usedAges[w]) < targetAges[w]*p[w] for w in wineClasses)
			for w in wineClasses
				for i in usedAges[w]
					prodSpecActions[w][i] += 1
				end
			end
			return p,prodSpecActions
		end
	end
	return nothing
end
#------------------------------------------------#
#checks whether a state-and-action pair is compatible
function isCompatible(state,action,totalDemand)
	if any(state[k] < action[k] for k in 1:length(action))
		return false
	elseif action[end] < min(totalDemand, state[end])
		return false
	# elseif sum(state)*(1/3) < sum(action)
	# 	return false
	else
		return true
	end
end
#------------------------------------------------#
function getInitTransitions(seed,newState,numAges,combinedDecay,transState=vcat(zeros(numAges-1),[1.0],zeros(numAges-1)))
	currState = copy(transState)
	currentElement = copy(seed)
	if currentElement >= numAges
		return currState'
	end
	transitionSpace = Array{Float64}(undef,0,numAges*2-1)
	prevInv = newState[currentElement]
	if prevInv > 0
		for d in 0:prevInv
			transState = copy(currState)
			transState[currentElement] = prevInv - d; transState[numAges+currentElement]=d; transState[numAges] = transState[numAges]*combinedDecay[prevInv][currentElement+1,d+1]
			newTransitions = getInitTransitions(currentElement+1,newState,numAges,combinedDecay,transState)
			transitionSpace = vcat(transitionSpace,newTransitions)
		end
	else
		transState = copy(currState)
		newTransitions = getInitTransitions(currentElement+1,newState,numAges,combinedDecay,transState)
		transitionSpace = vcat(transitionSpace,newTransitions)
	end
	return transitionSpace
end
#------------------------------------------------#
function getDecays(seed,newState,numAges,combinedDecay,stateMultipliers,supermarketContribution,transState=[1.0,1.0,0.0])
	currState = copy(transState)
	currentElement = copy(seed)
	if currentElement > numAges
		return currState'
	end
	transitionSpace = Array{Float64}(undef,0,3)
	prevInv = newState[currentElement]
	if prevInv > 0
		for d in 0:prevInv
			transState = copy(currState)
			transState[1] += (prevInv - d)*stateMultipliers[currentElement]; transState[2] = transState[2]*combinedDecay[prevInv][currentElement,d+1]; transState[3] += d*supermarketContribution[currentElement]
			newTransitions = getDecays(currentElement+1,newState,numAges,combinedDecay,stateMultipliers,supermarketContribution,transState)
			transitionSpace = vcat(transitionSpace,newTransitions)
		end
	else
		transState = copy(currState)
		newTransitions = getDecays(currentElement+1,newState,numAges,combinedDecay,stateMultipliers,supermarketContribution,transState)
		transitionSpace = vcat(transitionSpace,newTransitions)
	end
	return transitionSpace
end
#------------------------------------------------#
#checks whether a blending state is a possible transition state of a particular purchasing state
function isTransState(input,tState)
	for i in 1:length(input)
		if tState[i] > input[i]
			return false
		end
	end
	return true
end
#------------------------------------------------#
function getBlending(action)
	a = copy(action)
	occupiedPositions = 0
	blendingPattern = []
	for i in 1:length(a)
		while a[i] > 0
			occupiedPositions += 1
			append!(blendingPattern,i)
			a[i] -= 1
		end
	end
	return blendingPattern
end
#------------------------------------------------#
function remove!(a, item)
    deleteat!(a, findfirst(x->x==item, a))
end
#------------------------------------------------#
integratedBlendingAndPurchasing()
