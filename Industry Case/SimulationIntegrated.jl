using JSON
using Distributions
using Dates
using Base.Threads
using Random
using StatsBase
using DataFrames
using Statistics
using XLSX
#-------------------------------------------#

#-------------------------------------------#
function simulationIntegrated()
	for nP in [2]
		for sF in [25/6.0]
			for pU in [true]
				pUString = pU ? "PU" : "NoPU"
				pSrange = pU ? [20] : ["na"]
				for sprRev in [3]
					for oDP in [0.25]
						for pS in pSrange

							start = now()

							showRule = false

							ages = [i for i in 1:6]
							numAges = ages[end]

							scalingFactor = sF

							agesLarge = [i for i in 1:25]
							numAgesLarge = agesLarge[end]

							largeSpread = 10

							ageBuckets = [1,1,1,1,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6]

							#approximation case
							BAp = false
							SAp = false; protectionLevel = 6
							PrAp = false
							PuAp = false
							PuSAp = false
							PuPrAp = false
							TAp = false

							if PuAp || TAp || BAp || PrAp || SAp || PuSAp ||PuPrAp || showRule
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
								simulationInstance = simulationInstance * string("/Analysis/Simulation_Aggregated/",nP,"nP_",sF,"sF_",pUString,"pU_",sprRev,"sprRev_",oDP,"oDP_",pS,"pS_Simulation.json")
								open(simulationInstance, "r") do f
									pResults = JSON.parse(f)
								end
								purchasingResults = pResults["priceStockPurchasing"]
								purchasing3Results = pResults["priceStock3Purchasing"]
								for p in keys(purchasingResults) for k in keys(purchasingResults[p]) if purchasingResults[p][k] == nothing purchasingResults[p][k] = -1 end end end
								for p in keys(purchasing3Results) for k in keys(purchasing3Results[p]) if purchasing3Results[p][k] == nothing purchasing3Results[p][k] = -1 end end end
								optimalBlendingInv = pResults["ageStructureBlending"]
								optimalPurchasingInv = vcat([pResults["averagePurchasing"]],pResults["ageStructurePurchasing"])
								optimalPostBlendingInv = pResults["ageStructurePurchasing"]
							end

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
							if pU
								nPriceLevels = 9
							else
								nPriceLevels = 1
							end
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
								sort!(harvestCosts)
							else
								harvestCosts = [wineCostMean]
							end
							HtoPC = round(brandContribution[end]/mean(harvestCosts), digits=1)
							println(harvestCosts)

							#maximally allowed spread in blending
							maxSpread = 2

							#indicator for number of harvest price scenarios
							minYield = 0

							usePurchasingPredictions = false
							useBlendingPredictions = false

							start = now()

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
								for k in 1:currLength
									for i in 0:detDemand[p]
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
								yieldProbability = [untruncYield[i]/sum(values(untruncYield)) for i in 0:nPriceLevels-1]
								println(yieldProbability, sum(values(yieldProbability)))
							else
								yieldProbability = [1]
							end
							yieldItems = Array{Int,1}(collect(minYield:nPriceLevels-1))

							#load results from integrated model
							results = Dict()
							if pU
								currentInstance = string("",nP,"nP_",sF,"sF_","PU","pU_",sprRev,"sprRev_",oDP,"oDP_",pS,"pS")
							else
								currentInstance = string("",nP,"nP_",sF,"sF_","NoPU","pU_",sprRev,"sprRev_",oDP,"oDP_",pS,"pS")
							end
							open(string(string(@__DIR__),"/Results/",currentInstance,"_MDP.json"), "r") do f
								results = JSON.parse(f)
							end
							if TAp
								TApresults = Dict()
								open(string(string(@__DIR__),"/Results/",currentInstance,"_TAp.json"), "r") do f
									TApresults = JSON.parse(f)
								end
							end

							stateSpaceSize = length(values(results["stateSpace"]))
							purchasingStateSpaceSize = length(values(results["purchasingStateSpace"]))
							stateSpace = Dict(parse(Int32,i)=>results["stateSpace"][i] for i in keys(results["stateSpace"]))
							purchasingStateSpace = Dict(parse(Int32,i)=>results["purchasingStateSpace"][i] for i in keys(results["purchasingStateSpace"]))
							blendingStateAge = [sum(stateSpace[s][i]*i for i in 1:numAges)/sum(stateSpace[s]) for s in 1:stateSpaceSize]
							purchasingStateAge = [sum(purchasingStateSpace[s][i]*i for i in 2:numAges)/max(1,sum(purchasingStateSpace[s][2:end])) for s in 1:purchasingStateSpaceSize]
							inventoryLevelBlending = [sum(stateSpace[s]) for s in 1:stateSpaceSize]
							inventoryLevelPurchasing = [sum(purchasingStateSpace[s][2:end]) for s in 1:purchasingStateSpaceSize]
							purchasingSpace = [i for i in 0:nInventoryLevels-1]
							blendingSpace = Dict(parse(Int16,i)=>results["blendingActionSpace"][i] for i in keys(results["blendingActionSpace"]))
							blendingSpaceSize = length(values(blendingSpace))
							fulfillProduct = Dict(a => Dict(w => Int(results["fulfillProduct"][string(a)][string(w)]) for w in wineClasses) for a in 1:blendingSpaceSize)

							if !useBlendingPredictions
								if TAp
									optimalBlending = Dict(parse(Int32,s)=>TApresults["blending"][s] for s in keys(results["optimalBlending"]))
								else
									optimalBlending = Dict(parse(Int32,s)=>results["optimalBlending"][s] for s in keys(results["optimalBlending"]))
								end
							else
								res = Dict()
								open(string(string(@__DIR__),"/DecisionTree/PurchasingAmountTrees/",currentInstance,"_Predictions.json"), "r") do f
									res = JSON.parse(f)
								end
								optB = res["blendingPredictions"]
								optimalBlending = Dict(i => optB[i] for i in 1:stateSpaceSize)
							end
							if !usePurchasingPredictions
								if TAp
									optimalPurchasing = Dict(parse(Int32,s)=>TApresults["purchasing"][s] for s in keys(results["optimalPurchasing"]))
								else
									optimalPurchasing = Dict(parse(Int32,s)=>results["optimalPurchasing"][s] for s in keys(results["optimalPurchasing"]))
								end
							else
								res = Dict()
								open(string(string(@__DIR__),"/DecisionTree/PurchasingAmountTrees/",currentInstance,"_Predictions.json"), "r") do f
									res = JSON.parse(f)
								end
								optP = res["purchasingPredictions"]
								optimalPurchasing = Dict(i=>Int(round(optP[i],digits=0)) for i in 1:purchasingStateSpaceSize)
							end

							blendingProduct = Dict(a => Dict(w => results["blendingProduct"][string(a)][string(w)] for w in wineClasses) for a in 1:blendingSpaceSize)
							blendingProductRev = Dict(w => [blendingProduct[a][w] for a in 1:blendingSpaceSize] for w in wineClasses)
							prodSpecActions = Dict(a => Dict(w => results["prodSpecActions"][string(a)][string(w)] for w in wineClasses) for a in 1:blendingSpaceSize)

							averageReward = results["averageReward"]

							averageAge = Dict(a=>Dict(w=>((prodSpecActions[a][w])' * ages)/max(fulfillProduct[a][w],1) for w in wineClasses) for a in 1:blendingSpaceSize)
							stateBlendingActions = Dict(parse(Int32,s)=>results["stateBlendingActions"][s] for s in keys(results["stateBlendingActions"]))

							stateOutdating = [stateSpace[s][end] - blendingSpace[optimalBlending[s]][end] for s in 1:stateSpaceSize]
							blendingContribution = [sum(fulfillProduct[a][w]*(brandContribution[w] - holdingCosts*(averageAge[a][w]-targetAges[w])) for w in wineClasses) for a in 1:blendingSpaceSize]
							blendingStateContribution = [blendingContribution[optimalBlending[s]] for s in 1:stateSpaceSize]
							afterBlendingState = Dict(s=>Int(((stateSpace[s]-blendingSpace[optimalBlending[s]])[1:end-1])' * stateMultipliers[2:end] + 1) for s in 1:stateSpaceSize)
							reversedStateSpace = Dict(v=>k for (k,v) in stateSpace)
							stateHarvestCosts = [harvestCosts[purchasingStateSpace[s][1]+1] for s in 1:purchasingStateSpaceSize]
							statePurchasingCosts = [stateHarvestCosts[s] * optimalPurchasing[s] + holdingCosts*(optimalPurchasing[s]+sum(purchasingStateSpace[s][2:end])) for s in 1:purchasingStateSpaceSize]
							afterPurchasingState = Dict(s => 1 + optimalPurchasing[s]*stateMultipliers[1] + (stateMultipliers[2:end])' * purchasingStateSpace[s][2:end] for s in 1:purchasingStateSpaceSize)

							blendingPatterns = Dict(w => collect(Set(copy(blendingProductRev[w]))) for w in wineClasses)
							for w in wineClasses
								remove!(blendingPatterns[w], [])
							end

							noBlendingActions = Dict(w=>[b for b in 1:blendingSpaceSize if !isBlend(blendingProduct[b][w]) && fulfillProduct[b][w]>=2] for w in wineClasses)
							println("noBlendingActions: ",noBlendingActions)

							blendingPatternActions = Dict(w=>Dict(p=>[] for p in blendingPatterns[w]) for w in wineClasses)
							for b in 1:blendingSpaceSize
								for w in wineClasses
									if !isempty(blendingProduct[b][w])
										append!(blendingPatternActions[w][blendingProduct[b][w]], b)
									end
								end
							end

							maxStateSales = Dict(s=>Dict(a=>0 for a in ages) for s in 1:stateSpaceSize)
							stateBlendingProportion = Dict(s=>Dict(a=>0.0 for a in ages) for s in 1:stateSpaceSize)
							stateBlendingOptions = Dict(s => Dict(w => [] for w in wineClasses) for s in 1:stateSpaceSize)
							for s in 1:stateSpaceSize
								for b in stateBlendingActions[s]
									for w in wineClasses
										if !isempty(blendingProduct[b][w]) && !(blendingProduct[b][w] in stateBlendingOptions[s][w])
											append!(stateBlendingOptions[s][w],tuple(blendingProduct[b][w]))
										end
									end
									for a in ages
										if blendingSpace[b][a] > maxStateSales[s][a]
											maxStateSales[s][a] = blendingSpace[b][a]
										end
									end
								end
								for a in ages
									if maxStateSales[s][a] > 0
										stateBlendingProportion[s][a] = blendingSpace[optimalBlending[s]][a] / maxStateSales[s][a]
									else
										stateBlendingProportion[s][a] = -1.0
									end
								end
							end


							decayItems = Dict(); decayWeights = Dict()
							for a in 1:numAges
								decayItems[a]=Dict(); decayWeights[a]=Dict()
								for s in 1:nInventoryLevels-1
									decayItems[a][s] = Array{Int,1}([]); decayWeights[a][s]=Array{Float64,1}([])
									for d in 0:s
										append!(decayItems[a][s],d); append!(decayWeights[a][s],combinedDecay[s][a,d+1])
									end
								end
							end
							# println(combinedDecay)
							# fig,ax = subplots(1,1)
							# for k in 1:ages[end]
							# 	ax.plot(k,decayProbability[k],"dw", ms=8, label="discrete_weibull_pmf"*string(k))
							# 	ax.vlines(k,0,decayProbability[k],colors="b", lw=5, alpha=0.5)
							# end
							# show()
							final = now()-start
							println("Time to set up data and probabilities: ",final)

							#SIMULATION STUDY

							#set number of iterations
							numIterations = 20_000_000
							#find appropriate starting state and initialize data
							startingState = rand(1:stateSpaceSize)
							while sum(stateSpace[startingState]) < totalDemand*4
								startingState = rand(1:stateSpaceSize)
							end
							warm_up_period = 1_000_000
							warm_up_iteration = 0
							while warm_up_iteration <= warm_up_period
								warm_up_iteration += 1
								iterationData = performIteration(afterBlendingState[startingState],afterPurchasingState,yieldItems,Array{Float64,1}(yieldProbability),decayItems,decayWeights,numAges,stateSpace,supermarketContribution,stateMultipliers)
								startingState = iterationData[1]
							end

							if showRule
								#load decision trees
								for i in 1:10
									trees = Dict()
									treefile = string("/DecisionTree/PurchasingAmountTrees/",currentInstance,"_Trees.json")
									open(string(string(@__DIR__),treefile), "r") do f
										trees = JSON.parse(f)
									end
									purchasingTree = trees["6"]["true"]["purchasing"]
									blendingFFTree = Dict(c=> trees["6"]["true"]["blendingFF"][string(c)] for c in wineClasses)

									pFeatureList = trees["pFeatureList"]
									bFeatureList = trees["bFeatureList"]

									blendingDetTree = Dict(c => Dict() for c in wineClasses)
									for c in wineClasses
										for i in 1:detDemand[c]
											blendingDetTree[c][i] = trees["6"]["true"]["blendingDet"*string(i)*string(c)]
										end
									end
									performInterpretableIteration(afterBlendingState[startingState],afterPurchasingState,yieldItems,Array{Float64,1}(yieldProbability),decayItems,decayWeights,numAges,stateSpace,supermarketContribution,stateMultipliers,detDemand,purchasingTree,blendingFFTree,blendingDetTree,fulfillProduct,blendingContribution,optimalBlendingInv,optimalPurchasingInv,harvestCosts,stateBlendingActions,wineClasses)
									iterationData = performIteration(afterBlendingState[startingState],afterPurchasingState,yieldItems,Array{Float64,1}(yieldProbability),decayItems,decayWeights,numAges,stateSpace,supermarketContribution,stateMultipliers)
									startingState = iterationData[1]
								end
							end

							println("warm up completed!")
							println("starting state: ",stateSpace[startingState])
							println("first carry on after sales: ",afterBlendingState[startingState])
							iteration = 0
							stateVisitsBlending = zeros(stateSpaceSize)
							stateVisitsPurchasing = zeros(purchasingStateSpaceSize)
							blendingFrequency = zeros(blendingSpaceSize)
							purchasingFrequency = Dict(i=>0 for i in 0:nInventoryLevels-1)
							purchasingFrequencyPrice = Dict(p=>Dict(i=> 0 for i in 0:nInventoryLevels-1) for p in 0:nPriceLevels-1)
							blendingStatesSequence = []
							purchasingStatesSequence = []
							decayedItemsSequence = []
							decayProfitsSequence = []
							annualProfits = []
							while iteration < numIterations
								iteration += 1
								iterationData = performIteration(afterBlendingState[startingState],afterPurchasingState,yieldItems,Array{Float64,1}(yieldProbability),decayItems,decayWeights,numAges,stateSpace,supermarketContribution,stateMultipliers)
								stateVisitsBlending[iterationData[1]] += 1
								stateVisitsPurchasing[iterationData[2]] += 1
								blendingFrequency[optimalBlending[iterationData[1]]] += 1
								purchasingFrequency[optimalPurchasing[iterationData[2]]] += 1
								purchasingFrequencyPrice[purchasingStateSpace[iterationData[2]][1]][optimalPurchasing[iterationData[2]]] += 1
								append!(blendingStatesSequence,iterationData[1])
								append!(purchasingStatesSequence,iterationData[2])
								append!(decayProfitsSequence,iterationData[3])
								append!(decayedItemsSequence,iterationData[4])
								append!(annualProfits, -statePurchasingCosts[iterationData[2]] + iterationData[3] + blendingStateContribution[iterationData[1]])
								startingState = iterationData[1]
								if iteration%1_000_000 == 0
									println("", iteration, " iterations done")
									println("average profits: ", mean(annualProfits))
								end
							end

							# nShortTermRuns = 50
							# println("Start short-term runs from the following state: ", stateSpace[startingState])
							# shortTermStart = startingState
							# shortTermLength = 30
							# nRuns = 0
							# shortTermProfits = Dict(i=>[] for i in 1:nShortTermRuns)
							# while nRuns < nShortTermRuns
							# 	nRuns += 1
							# 	iteration = 0
							# 	startingState = shortTermStart
							# 	while iteration < shortTermLength
							# 		iteration += 1
							# 		iterationData = performIteration(afterBlendingState[startingState],afterPurchasingState,yieldItems,Array{Float64,1}(yieldProbability),decayItems,decayWeights,numAges,stateSpace,supermarketContribution,stateMultipliers)
							# 		append!(shortTermProfits[nRuns], -statePurchasingCosts[iterationData[2]] + iterationData[3] + blendingStateContribution[iterationData[1]])
							# 		startingState = iterationData[1]
							# 	end
							# end
							# println("Short-term average profits: ")
							# for n in 1:nShortTermRuns
							# 	println(mean(shortTermProfits[n]))
							# end
							#
							# nShortTermRuns = 50
							# badState = reversedStateSpace[[0,0,0,2,0,5,5]]
							# println("Start short term runs from bad blending state: ",stateSpace[badState]," with optimal policy ", optimalBlending[badState])
							# nRuns = 0
							# shortTermProfits = Dict(i=>[] for i in 1:nShortTermRuns)
							# while nRuns < nShortTermRuns
							# 	nRuns += 1
							# 	iteration = 0
							# 	startingState = badState
							# 	while iteration < shortTermLength
							# 		iteration += 1
							# 		iterationData = performIteration(afterBlendingState[startingState],afterPurchasingState,yieldItems,Array{Float64,1}(yieldProbability),decayItems,decayWeights,numAges,stateSpace,supermarketContribution,stateMultipliers)
							# 		append!(shortTermProfits[nRuns], -statePurchasingCosts[iterationData[2]] + iterationData[3] + blendingStateContribution[iterationData[1]])
							# 		startingState = iterationData[1]
							# 	end
							# end
							# println("Short-term average profits: ")
							# for n in 1:nShortTermRuns
							# 	println(mean(shortTermProfits[n]))
							# end

							#evaluate results
							println(purchasingFrequency)
							non_visits_blending = 0
							non_visits_purchasing = 0
							deliberateUnderfulfillment = Dict(c=>[] for c in wineClasses)
							for s in 1:stateSpaceSize
								if stateVisitsBlending[s] == 0
									non_visits_blending += 1
								end
								if s <= purchasingStateSpaceSize && stateVisitsPurchasing[s] == 0
									non_visits_purchasing += 1
								end
								for c in wineClasses
									if fulfillProduct[optimalBlending[s]][c] < maximum(fulfillProduct[f][c] for f in stateBlendingActions[s])
										append!(deliberateUnderfulfillment[c],s)
									end
								end
							end
							println(averageReward)
							stockLevel = zeros((numAges-1)*(nInventoryLevels-1) + 1); stockLevelPrice3 = zeros(nPriceLevels,(nInventoryLevels-1)*2 + 1); stockLevelPrice = zeros(nPriceLevels,(numAges-1)*nInventoryLevels + 1); stockLevelCount = zeros((numAges)*(nInventoryLevels-1) + 1); for s in 1:purchasingStateSpaceSize stockLevel[inventoryLevelPurchasing[s]+1] += stateVisitsPurchasing[s]; stockLevelPrice3[purchasingStateSpace[s][1]+1,sum(purchasingStateSpace[s][2:3])+1] += stateVisitsPurchasing[s]; stockLevelPrice[purchasingStateSpace[s][1]+1,inventoryLevelPurchasing[s]+1] += stateVisitsPurchasing[s]; stockLevelCount[inventoryLevelPurchasing[s]+1] += 1 end

							orderUpTo = Dict(p => zeros(numAges*nInventoryLevels + 1) for p in 0:nPriceLevels-1)
							orderUpTo3 = Dict(p => zeros(4*nInventoryLevels + 1) for p in 0:nPriceLevels-1)

							stockLevelBlending =zeros(numAges*(nInventoryLevels-1) + 1); stockLevelBlendingCount = zeros(numAges*(nInventoryLevels-1) + 1); for s in 1:stateSpaceSize stockLevelBlending[inventoryLevelBlending[s]+1] += stateVisitsBlending[s]; stockLevelBlendingCount[inventoryLevelBlending[s]+1] += 1 end
							afterPurchasingInventory = zeros(numAges*(nInventoryLevels-1)+1); for s in 1:purchasingStateSpaceSize afterPurchasingInventory[inventoryLevelBlending[afterPurchasingState[s]]+1] += stateVisitsPurchasing[s] end
							afterPurchasingInventory3 = zeros(4*(nInventoryLevels-1)+1); for s in 1:purchasingStateSpaceSize afterPurchasingInventory3[sum(stateSpace[afterPurchasingState[s]][1:4])+1] += stateVisitsPurchasing[s] end
							blendingInventory3_end = zeros((numAges-2)*(nInventoryLevels-1)+1); for s in 1:stateSpaceSize blendingInventory3_end[sum(stateSpace[s][3:end])+1] += stateVisitsBlending[s] end
							blendingInventory4_end = zeros((numAges-3)*(nInventoryLevels-1)+1); for s in 1:stateSpaceSize blendingInventory4_end[sum(stateSpace[s][4:end])+1] += stateVisitsBlending[s] end
							blendingInventory5_end = zeros((numAges-4)*(nInventoryLevels-1)+1); for s in 1:stateSpaceSize blendingInventory5_end[sum(stateSpace[s][5:end])+1] += stateVisitsBlending[s] end

							println("stock levels purchasing: ", stockLevel, " count: ", stockLevelCount)
							println("stock levels blending: ", stockLevelBlending, " count: ", stockLevelBlendingCount)
							stockLevelsInvestigated = [i for i in 0:length(stockLevel)-1]
							avgAgeIntervals = [round(1+0.1*i, digits=1) for i in 0:(numAges-1)*10]
							avgAgeFrequencyPurchasing = zeros(length(avgAgeIntervals))
							for s in 1:purchasingStateSpaceSize
								if purchasingStateAge[s] >= 1
									c =Int(round(purchasingStateAge[s],digits=1)*10-9)
									avgAgeFrequencyPurchasing[c] += stateVisitsPurchasing[s]
								end
							end
							for i in 1:length(avgAgeIntervals)
								println(i, " avg age frequency: ", avgAgeFrequencyPurchasing[i])
							end
							purchasingStockPrice = Dict(p=>Dict(l=>[] for l in stockLevelsInvestigated) for p in 0:nPriceLevels-1)
							purchasingStock3Price = Dict(p=>Dict(l=>[] for l in 0:(nInventoryLevels-1)*3) for p in 0:nPriceLevels-1)

							for s in 1:purchasingStateSpaceSize
								ps = purchasingStateSpace[s]
								orderUpTo[ps[1]][inventoryLevelPurchasing[s]+optimalPurchasing[s]] += stateVisitsPurchasing[s]
								orderUpTo3[ps[1]][sum(ps[2:4])+optimalPurchasing[s]+1] += stateVisitsPurchasing[s]
								if inventoryLevelPurchasing[s] in stockLevelsInvestigated
									append!(purchasingStockPrice[ps[1]][inventoryLevelPurchasing[s]], optimalPurchasing[s]*stateVisitsPurchasing[s])
									append!(purchasingStock3Price[ps[1]][sum(ps[2:4])], optimalPurchasing[s]*stateVisitsPurchasing[s])
								end
							end
							meanPurchasingStockPrice = Array{Float64,2}(undef, nPriceLevels, length(stockLevelsInvestigated))
							meanPurchasingStock3Price = Array{Float64,2}(undef, nPriceLevels, (nInventoryLevels-1)*2+1)
							fill!(meanPurchasingStockPrice, 0.0); fill!(meanPurchasingStock3Price, 0.0);
							for p in 0:nPriceLevels-1
								orderUpTo[p] = orderUpTo[p]/sum(orderUpTo[p])
								orderUpTo3[p] = orderUpTo3[p]/sum(orderUpTo3[p])
								for l in stockLevelsInvestigated
									if purchasingStockPrice[p][l] != [] && stockLevelPrice[p+1,l+1] != 0
										meanPurchasingStockPrice[p+1,l+1]= sum(purchasingStockPrice[p][l])/stockLevelPrice[p+1,l+1]
									else
										meanPurchasingStockPrice[p+1,l+1]=-1
									end
								end
								for l in 0:(nInventoryLevels-1)*2
									if purchasingStock3Price[p][l] != [] && stockLevelPrice3[p+1,l+1] != 0
										meanPurchasingStock3Price[p+1,l+1]= sum(purchasingStock3Price[p][l])/stockLevelPrice3[p+1,l+1]
									else
										meanPurchasingStock3Price[p+1,l+1]=-1
									end
								end
							end

							timesAgeSold = zeros(numAges); for b in 1:blendingSpaceSize for i in ages timesAgeSold[i] += blendingFrequency[b]*blendingSpace[b][i] end end
							timesAgesSoldProduct = Dict(w => zeros(numAges) for w in wineClasses); for b in 1:blendingSpaceSize for i in ages for w in wineClasses timesAgesSoldProduct[w][i] += blendingFrequency[b]*prodSpecActions[b][w][i] end end end
							averageSelling = [timesAgeSold[i]/numIterations for i in 1:numAges]
							averageSellingProduct = Dict(w=> [timesAgesSoldProduct[w][i]/numIterations for i in 1:numAges] for w in wineClasses)
							avgTargetAgeExcess = Float64(averageSelling' * ages)/(sum(sum(averageSellingProduct[w])*targetAges[w] for w in wineClasses)) 
							avgTargetAgeExcessProduct = Dict(w=>((averageSellingProduct[w])' * ages)/(sum(averageSellingProduct[w])*targetAges[w]) for w in wineClasses)
							
							onTargetWide = Dict(w=> [p for p in blendingPatterns[w] if mean(p) == targetAges[w] && p[1] != p[end] && sum(p) == detDemand[w]*targetAges[w]] for w in wineClasses)
							println("on target wide: ", onTargetWide)
							onTargetTight = Dict(w=> [p for p in blendingPatterns[w] if mean(p) == targetAges[w] && p[1] == p[end] && sum(p) == detDemand[w]*targetAges[w]] for w in wineClasses)
							println("on target tight: ", onTargetTight)

							preferTight = Dict(w => sum([stateVisitsBlending[s] for s in 1:stateSpaceSize if optimalBlending[s] in vcat([blendingPatternActions[w][p] for p in onTargetTight[w]]...) && any(a in stateBlendingActions[s] for a in vcat([blendingPatternActions[w][p] for p in onTargetWide[w]]...))]) for w in wineClasses)
							preferWide = Dict(w => sum([stateVisitsBlending[s] for s in 1:stateSpaceSize if optimalBlending[s] in vcat([blendingPatternActions[w][p] for p in onTargetWide[w]]...) && any(a in stateBlendingActions[s] for a in vcat([blendingPatternActions[w][p] for p in onTargetTight[w]]...))]) for w in wineClasses)
							println(preferTight)
							println(preferWide)
							percentageWide = Dict(w => preferWide[w]/(preferWide[w]+preferTight[w]) for w in wineClasses)
							println(percentageWide)

							println("Relative use of age in sales: ", string(timesAgeSold/(numIterations*totalDemand)))

							mostVisitsBlending = partialsortperm(stateVisitsBlending,1:10,rev=true)
							mostVisitsPurchasing = partialsortperm(stateVisitsPurchasing,1:10,rev=true)
							mostBlendingActions = partialsortperm(blendingFrequency,1:10,rev=true)
							println("Most visited states blending: "); for s in mostVisitsBlending println(stateSpace[s], " ", stateVisitsBlending[s]) end
							println("Most visited states purchasing: "); for s in mostVisitsPurchasing println(stateSpace[s], " ", stateVisitsPurchasing[s]) end
							println("Most frequent blending actions: "); for b in mostBlendingActions println(blendingSpace[b], " ", blendingFrequency[b]) end

							blendingPatternsFrequency = Dict(c=>Dict(p=>sum(blendingFrequency[b] for b in blendingPatternActions[c][p]) for p in blendingPatterns[c]) for c in wineClasses)
							blendingAvailable = Dict(c=>Dict(p=>sum([stateVisitsBlending[s] for s in 1:stateSpaceSize if p in stateBlendingOptions[s][c]]) for p in blendingPatterns[c]) for c in wineClasses)
							println("blendingAvailable: ", blendingAvailable)

							avgDemandFF = Dict(w=>sum(fulfillProduct[optimalBlending[s]][w]*stateVisitsBlending[s] for s in 1:stateSpaceSize)/numIterations for w in wineClasses)
							averagePurchasing = sum(purchasingFrequency[i]*i for i in 0:nInventoryLevels-1)/numIterations
							avgBlendingInventory = [sum(stateSpace[s][i] * stateVisitsBlending[s] for s in 1:stateSpaceSize) for i in 1:numAges]/numIterations
							avgPurchasingInventory = [sum(purchasingStateSpace[s][i] * stateVisitsPurchasing[s] for s in 1:purchasingStateSpaceSize) for i in 2:numAges]/numIterations
							countAvailable = zeros(numAges)
							countProportion = zeros(numAges)
							for s in 1:stateSpaceSize
								for i in 1:numAges
									if stateBlendingProportion[s][i] >= 0
										countAvailable[i] += stateVisitsBlending[s]
										countProportion[i] += stateVisitsBlending[s]*stateBlendingProportion[s][i]
									end
								end
							end
							ageSellingIfAvailable = countProportion/countAvailable
							demandPatternFrequency = Dict(p=>0 for p in demandPatterns)

							for b in 1:blendingSpaceSize
								demandPatternFrequency[[fulfillProduct[b][w] for w in wineClasses]] += blendingFrequency[b]
							end
							println("non-visits-blending: ", non_visits_blending, " non-visits-purchasing: ", non_visits_purchasing)
							println("Demand Fulfillment: ", avgDemandFF)
							println("Demand Fulfillment Frequency: ", demandPatternFrequency)
							println("Average inventory level: in blending state: ", sum(avgBlendingInventory), " in purchasing state: ", sum(avgPurchasingInventory))
							println("Inventory age structure: blending: ", avgBlendingInventory, " purchasing: ", avgPurchasingInventory)
							println("Average reward comparison: MDP model: ", averageReward, ", Simulation study: ", mean(annualProfits))
							println("Standard deviation of rewards: ", std(annualProfits))
							println("Average selling per age: ", ageSellingIfAvailable)
							println("Frequency of blending patterns: ", blendingPatternsFrequency)
							println("Percentage of applied blending, when available: "); for c in wineClasses for p in blendingPatterns[c] println(c, " ", p, " ", blendingPatternsFrequency[c][p]/blendingAvailable[c][p]) end end
							fullDemandFF = Dict(w=>sum(blendingFrequency[b] for b in 1:blendingSpaceSize if fulfillProduct[b][w] == detDemand[w]) for w in wineClasses)
							println("Total number of full demand fulfillment: ",fullDemandFF)
							# println("No blending when available: Young: ", string(sum(blendingFrequency[b] for b in noBlendingActions["young"])/sum(stateVisitsBlending[s] for s in 1:stateSpaceSize if any(b in noBlendingActions["young"] for b in stateBlendingActions[s]))), " Old: ", string(sum(blendingFrequency[b] for b in noBlendingActions["old"])/sum(stateVisitsBlending[s] for s in 1:stateSpaceSize if any(b in noBlendingActions["old"] for b in stateBlendingActions[s]))))
							println("No blending frequency: ",Dict(w=>sum(blendingFrequency[b] for b in noBlendingActions[w]) for w in wineClasses))

							#println("Deliberate Underfulfillment qty: ", Dict(w=>sum(stateVisitsBlending[s] for s in deliberateUnderfulfillment[w]) for w in wineClasses))
							#println("Deliberate demand underfulfillment: ", Dict(w=>sum(stateVisitsBlending[s] for s in deliberateUnderfulfillment[w])/(numIterations-fullDemandFF[w]) for w in wineClasses))

							# println([stateSpace[s] for s in deliberateUnderfulfillment["old"]])
							pricePurchasing = Dict(p=>[purchasingFrequencyPrice[p][i]/sum(values(purchasingFrequencyPrice[p])) for i in 0:nInventoryLevels-1] for p in 0:nPriceLevels-1)
							priceStockPurchasing = Dict(p=>Dict(l=>meanPurchasingStockPrice[p+1,l+1] for l in stockLevelsInvestigated) for p in 0:nPriceLevels-1)
							priceStock3Purchasing = Dict(p=>Dict(l=>meanPurchasingStock3Price[p+1,l+1] for l in 0:(nInventoryLevels-1)*2) for p in 0:nPriceLevels-1)
							jsonRes = Dict()
							jsonRes["visitedBlending"] = stateVisitsBlending
							jsonRes["averageReward"] = mean(annualProfits)
							jsonRes["rewardSequence"] = annualProfits
							jsonRes["visitedPurchasing"] = stateVisitsPurchasing
							jsonRes["ageStructureBlending"] = avgBlendingInventory
							jsonRes["ageStructurePurchasing"] = append!(avgPurchasingInventory,avgBlendingInventory[end]-averageSelling[end])
							jsonRes["averagePurchasing"] = averagePurchasing
							jsonRes["averageSelling"] = averageSelling
							jsonRes["averageSellingProduct"] = averageSellingProduct
							jsonRes["avgTargetAgeExcess"] = avgTargetAgeExcess
							jsonRes["avgTargetAgeExcessProduct"] = avgTargetAgeExcessProduct
							jsonRes["percentageWide"] = percentageWide
							jsonRes["pricePurchasing"] = pricePurchasing
							jsonRes["serviceLevels"] = [avgDemandFF[w]/detDemand[w] for w in wineClasses]
							jsonRes["priceStockPurchasing"] = priceStockPurchasing
							jsonRes["priceStock3Purchasing"] = priceStock3Purchasing
							jsonRes["orderUpTo"] = orderUpTo
							jsonRes["orderUpTo3"] = orderUpTo3
							if usePurchasingPredictions == false
								xlsxPath = string(string(@__DIR__),"/Analysis/",currentInstance,"_Simulation.xlsx")
								jsonPath = string(string(@__DIR__),"/Analysis/",currentInstance,"_Simulation.json")
								open(jsonPath, "w") do f
									JSON.print(f,jsonRes,4)
								end
							else
								jsonPath = string(string(@__DIR__),"/Analysis/",currentInstance,"_Simulation_Prediction.json")
								xlsxPath = string(string(@__DIR__),"/Analysis/",currentInstance,"_Simulation_Prediction.xlsx")
								open(jsonPath, "w") do f
									JSON.print(f,jsonRes,4)
								end
							end
							XLSX.openxlsx(xlsxPath, mode="w") do xf
								OVS = xf[1]
								XLSX.rename!(OVS, "Overview")

								OVS["A1"] = "Average annual profit"
								OVS["B1"] = mean(annualProfits)
								OVS["A2"] = "Profit standard deviation"
								OVS["B2"] = std(annualProfits)

								OVS["A3"] = "Stock levels purchasing"
								OVS["B3"] = [i for i in 0:(numAges-1)*(nInventoryLevels-1)]
								OVS["A4"] = "Iteration count"
								OVS["B4"] = stockLevel
								OVS["A5"] = "State count"
								OVS["B5"] = stockLevelCount

								OVS["A7"] = "Stock levels blending"
								OVS["B7"] = [i for i in 0:(numAges)*(nInventoryLevels-1)]
								OVS["A8"] = "Iteration count"
								OVS["B8"] = stockLevelBlending
								OVS["A9"] = "State count"
								OVS["B9"] = stockLevelBlendingCount

								OVS["A11"] = "Average age intervals"
								OVS["B11"] = avgAgeIntervals
								OVS["A12"] = "Frequency"
								OVS["B12"] = avgAgeFrequencyPurchasing

								OVS["A14"] = "Age classes"
								OVS["B14"] = ages
								OVS["A15"] = "Relative use of age in sales: "
								OVS["B15"] = timesAgeSold/sum(timesAgeSold)
								OVS["A16"] = "Average selling if available"
								OVS["B16"] = [mean(ageSellingIfAvailable[a]) for a in 1:numAges]
								OVS["A17"] = "Average inventory level blending"
								OVS["B17"] = avgBlendingInventory
								OVS["A18"] = "Average inventory level purchasing"
								OVS["B18"] = avgPurchasingInventory
								OVS["A19"] = "Average sales per age class"
								OVS["B19"] = averageSelling
								OVS["A20"] = "Average purchasing amount"
								OVS["B20"] = averagePurchasing

								OVS["A21"] = "Stock levels after purchasing"
								OVS["B21"] = [i for i in 0:(numAges)*(nInventoryLevels-1)]
								OVS["A22"] = "Iteration count full inventory"
								OVS["B22"] = afterPurchasingInventory
								OVS["A23"] = "Iteration count age classes 1-3"
								OVS["B23"] = afterPurchasingInventory3
								OVS["A24"] = "Iteration count age classes 3-end"
								OVS["B24"] = blendingInventory3_end
								OVS["A25"] = "Iteration count age classes 4-end"
								OVS["B25"] = blendingInventory4_end

								OVS["A27"] = "Demand patterns"
								OVS["B27"] = [string(i) for i in demandPatterns]
								OVS["A28"] = "Frequency"
								OVS["B28"] = [demandPatternFrequency[p] for p in demandPatterns]

								OVS["A30"] = "Blending patterns 2yo"
								OVS["B30"] = [string(p) for p in blendingPatterns[1]]
								OVS["A31"] = "Frequency"
								OVS["B31"] = [sum(blendingFrequency[b] for b in blendingPatternActions[1][p]) for p in blendingPatterns[1]]
								OVS["A32"] = "Usage if available"
								OVS["B32"] = [blendingPatternsFrequency[1][p]/blendingAvailable[1][p] for p in blendingPatterns[1]]

								if nProducts > 2
									OVS["A33"] = "Blending patterns 3yo"
									OVS["B33"] = [string(p) for p in blendingPatterns[2]]
									OVS["A34"] = "Frequency"
									OVS["B34"] = [sum(blendingFrequency[b] for b in blendingPatternActions[2][p]) for p in blendingPatterns[2]]
									OVS["A34"] = "Usage if available"
									OVS["B34"] = [blendingPatternsFrequency[2][p]/blendingAvailable[2][p] for p in blendingPatterns[2]]
								end

								OVS["A36"] = "Blending patterns 4yo"
								OVS["B36"] = [string(p) for p in blendingPatterns[nProducts]]
								OVS["A37"] = "Frequency"
								OVS["B37"] = [sum(blendingFrequency[b] for b in blendingPatternActions[nProducts][p]) for p in blendingPatterns[nProducts]]
								OVS["A38"] = "Usage if available"
								OVS["B38"] = [blendingPatternsFrequency[nProducts][p]/blendingAvailable[nProducts][p] for p in blendingPatterns[nProducts]]

								OVS["A40"] = "Purchasing decisions"
								OVS["B40"] = [i for i in 0:nInventoryLevels-1]
								OVS["A41"] = "Frequency"
								OVS["B41"] = [purchasingFrequency[i] for i in 0:nInventoryLevels-1]
								OVS["A42"] = "Percentage"
								OVS["B42"] = [purchasingFrequency[i] for i in 0:nInventoryLevels-1]/sum(values(purchasingFrequency))

								OVS["A43"] = "Average outdating"
								OVS["B43"] = avgBlendingInventory[end] - averageSelling[end]

								# OVS["A44"] = "Average no. of decayed items"
								# OVS["B44"] = mean(decayedItemsSequence)

								OVS["A45"] = "Non-visited states blending"
								OVS["B45"] = non_visits_blending
								OVS["A46"] = "Non-visited states purchasing"
								OVS["B46"] = non_visits_purchasing

								OVS["A47"] = "Average demand fulfillment 2yo"
								OVS["B47"] = avgDemandFF[1]
								OVS["C47"] = avgDemandFF[1]/detDemand[1]

								if nProducts > 2
									OVS["A48"] = "Average demand fulfillment 3yo"
									OVS["B48"] = avgDemandFF[2]
									OVS["C48"] = avgDemandFF[2]/detDemand[2]
								end

								OVS["A49"] = "Average demand fulfillment 4yo"
								OVS["B49"] = avgDemandFF[nProducts]
								OVS["C49"] = avgDemandFF[nProducts]/detDemand[nProducts]

								if !usePurchasingPredictions
									OVS["A50"] = "Deliberate Underfulfillment"
									OVS["B50"] = [sum(stateVisitsBlending[s] for s in deliberateUnderfulfillment[c]) for c in wineClasses]
									OVS["B51"] = [sum(stateVisitsBlending[s] for s in deliberateUnderfulfillment[c])/(numIterations-fullDemandFF[c]) for c in wineClasses]
								end

								# OVS["A51"] = "No blending young"
								# OVS["B51"] = sum(blendingFrequency[b] for b in noBlendingActions["young"])
								# OVS["A52"] = "Percentage if available"
								# OVS["B52"] = sum(blendingFrequency[b] for b in noBlendingActions["young"])/sum(stateVisitsBlending[s] for s in 1:stateSpaceSize if any(b in noBlendingActions["young"] for b in stateBlendingActions[s]))
								# OVS["A53"] =  "No blending old"
								# OVS["B53"] = sum(blendingFrequency[b] for b in noBlendingActions["old"])
								# OVS["A54"] = "Percentage if available"
								# OVS["B54"] = sum(blendingFrequency[b] for b in noBlendingActions["old"])/sum(stateVisitsBlending[s] for s in 1:stateSpaceSize if any(b in noBlendingActions["old"] for b in stateBlendingActions[s]))

								OVS["A56"]="Most visited blending states"
								OVS["A57"]="Frequency"
								OVS["B56"]=[string(Array{Int64,1}(stateSpace[s])) for s in mostVisitsBlending]
								OVS["B57"]=[stateVisitsBlending[s] for s in mostVisitsBlending]

								OVS["A58"]="Most visited purchasing states"
								OVS["A59"]="Frequency"
								OVS["B58"]=[string(Array{Int64,1}(stateSpace[s])) for s in mostVisitsPurchasing]
								OVS["B59"]=[stateVisitsPurchasing[s] for s in mostVisitsPurchasing]

								OVS["A60"]="Most frequent blending actions"
								OVS["A61"]="Frequency"
								OVS["B60"]=[string(Array{Int64,1}(blendingSpace[b])) for b in mostBlendingActions]
								OVS["B61"]=[blendingFrequency[b] for b in mostBlendingActions]

								# OVS["A100"] = "Purchasing price and stock level based purchasing"
								# OVS["A101", dim=1] = harvestCosts
								# OVS["B100"] = stockLevelsInvestigated
								# if nProducts > 2
								# 	lastCol = "AX"
								# else
								# 	lastCol = "AX"
								# end
								# if pU
								# 	lastRow = "109"
								# else
								# 	lastRow = "101"
								# end
								# OVS[string("B101:",lastCol, lastRow)] = meanPurchasingStockPrice
								# #
								# if nProducts > 2
								# 	lastCol = "Z"
								# else
								# 	lastCol = "Z"
								# end
								# if pU
								# 	lastRow = "120"
								# else
								# 	lastRow = "112"
								# end
								# OVS["A111"] = "Purchasing price and stock level based purchasing"
								# OVS["A112", dim=1] = harvestCosts
								# OVS["B111"] = stockLevelsInvestigated
								# OVS[string("B112:",lastCol, lastRow)] = meanPurchasingStock3Price
								# OVS["B108"] = stockLevelYoung/youngTotal
								# OVS["A109"] = "Purchasing price and stock level based purchasing if avg age <= 3.5"
								# OVS["A110", dim=1] = harvestCosts
								# OVS["B109"] = stockLevelsInvestigated
								# OVS["B110:AF116"] = meanPurchasingStockPriceAge["young"]
								#
								# OVS["B117"] = stockLevelOld/oldTotal
								# OVS["A118"] = "Purchasing price and stock level based purchasing if avg age > 3.5"
								# OVS["A119", dim=1] = harvestCosts
								# OVS["B118"] = stockLevelsInvestigated
								# OVS["B119:AF125"] = meanPurchasingStockPriceAge["old"]

							end
						end
					end
			 	end
		  	end
	   	end
   end
end
#-------------------------------------------#
#perform simulation of one entire year (incl. blending/sales & purchasing)
function performIteration(afterBlendingState,afterPurchasingState,yieldItems,yieldProbability,decayItems,decayWeights,numAges,stateSpace,supermarketContribution,stateMultipliers,transState = 1)
	simulationState = copy(afterBlendingState)

	#simulate harvest price
	harvestPriceCategory = StatsBase.sample(yieldItems, Weights(yieldProbability))
	prePurchaseState = Int(simulationState + harvestPriceCategory*stateMultipliers[1])
	afterPurchaseState = afterPurchasingState[prePurchaseState]

	#simulate deccays
	decayProfits = 0; decayedItems = 0;
	for carryOnPos in 1:numAges
		currInv = stateSpace[afterPurchaseState][carryOnPos]
		if currInv > 0
			decay = StatsBase.sample(decayItems[carryOnPos][currInv], Weights(decayWeights[carryOnPos][currInv]))
			transState += (currInv - decay)*stateMultipliers[carryOnPos]; decayProfits += supermarketContribution[carryOnPos]*decay; decayedItems += decay
		end
	end
	return transState, prePurchaseState, decayProfits, decayedItems
end
#------------------------------------------------#
#perform simulation of one entire year (incl. blending/sales & purchasing)
function performInterpretableIteration(afterBlendingState,afterPurchasingState,yieldItems,yieldProbability,decayItems,decayWeights,numAges,stateSpace,supermarketContribution,stateMultipliers,detDemand,purchasingTree,blendingFFTree,blendingDetTree,fulfillProduct,blendingContribution,optimalBlendingInv,optimalPurchasingInv,harvestCosts,stateBlendingActions,wineClasses)
	simulationState = stateSpace[copy(afterBlendingState)]

	ages = [i for i in 1:numAges]
	optBAge = (optimalBlendingInv' * ages)/sum(optimalBlendingInv)
	optPAge = ((optimalPurchasingInv[2:5])' * (ages[2:end]))/sum(optimalPurchasingInv[2:end])

	features = Dict()
	optPFeature = Dict()

	features["pInv"] = sum(simulationState[2:5]-optimalPurchasingInv[2:5])/(numAges-1); optPFeature["pInv"] = sum(optimalPurchasingInv[2:5])/(numAges-1)
	for i in 2:numAges
		features[string("pInv",i)] = simulationState[i]-optimalPurchasingInv[i]; optPFeature[string("pInv",i)] = optimalPurchasingInv[i]
	end
	features["pInv23"] = sum(simulationState[2:3]-optimalPurchasingInv[2:3])/(2); optPFeature["pInv23"] = sum(optimalPurchasingInv[2:3])/(2)
	features["pInv234"] = sum(simulationState[2:4]-optimalPurchasingInv[2:4])/(3); optPFeature["pInv234"] = sum(optimalPurchasingInv[2:4])/(3)
	features["pInv34"] = sum(simulationState[3:4]-optimalPurchasingInv[3:4])/(2); optPFeature["pInv34"] = sum(optimalPurchasingInv[3:4])/(2)
	features["pInv345"] = sum(simulationState[3:5]-optimalPurchasingInv[3:5])/(3); optPFeature["pInv345"] = sum(optimalPurchasingInv[3:5])/(3)
	features["pInv45"] = sum(simulationState[4:5]-optimalPurchasingInv[4:5])/(2); optPFeature["pInv45"] = sum(optimalPurchasingInv[4:5])/(2)
	for c in wineClasses
		features[string("pFFMax",c)] = maximum(fulfillProduct[a][c] for a in stateBlendingActions[afterBlendingState])/detDemand[c]
	end
	features["pMaxRev"] = maximum(blendingContribution[a] for a in stateBlendingActions[afterBlendingState])
	features["pAge"] = (((simulationState[2:end])' * (ages[2:end]))/max(1,sum(simulationState[2:end]))-optPAge)/optPAge

	#simulate harvest price
	harvestPriceCategory = StatsBase.sample(yieldItems, Weights(yieldProbability))
	prePurchaseState = Int(afterBlendingState + harvestPriceCategory*stateMultipliers[1])

	features["pPrice"] = harvestCosts[harvestPriceCategory]

	println("Purchasing State: ", stateSpace[prePurchaseState])
	println("Optimal Purchasing Inventory: ", optimalPurchasingInv)

	pPrediction = nothing
	foundLeaf = false
	currentNode = "1"
	counter = 0
	while !foundLeaf
		counter += 1
		println(counter)
		if purchasingTree[currentNode]["isLeaf"]
			pPrediction = Int(round(purchasingTree[currentNode]["class"],digits=0))
			foundLeaf=true
		else
			currentFeature = purchasingTree[currentNode]["feature"]
			println("Feature: ", currentFeature)
			value = features[currentFeature]
			println("Value of current state: ", value)
			if currentFeature in ["pInv", "pInv2", "pInv3", "pInv4", "pInv5", "pInv23", "pInv34", "pInv45", "pInv234", "pInv345"]
				treeValue = purchasingTree[currentNode]["value"] - optPFeature[currentFeature]
			else
				treeValue = purchasingTree[currentNode]["value"]
			end
			println("Value of tree feature: ", treeValue)
			currentNode = value < treeValue ? string(purchasingTree[currentNode]["children"][1]) : string(purchasingTree[currentNode]["children"][2])
		end
	end
	println("Rule: ", pPrediction)
	afterPurchaseState = Int(afterBlendingState + pPrediction*stateMultipliers[1])

	#simulate deccays
	transState = 1
	decayProfits = 0; decayedItems = 0;
	for carryOnPos in 1:numAges
		currInv = stateSpace[afterPurchaseState][carryOnPos]
		if currInv > 0
			decay = StatsBase.sample(decayItems[carryOnPos][currInv], Weights(decayWeights[carryOnPos][currInv]))
			transState += (currInv - decay)*stateMultipliers[carryOnPos]; decayProfits += supermarketContribution[carryOnPos]*decay; decayedItems += decay
		end
	end

	simulationState = stateSpace[transState]
	println("Blending State: ", simulationState)

	optBFeature = Dict()

	features["bInv"] = sum(simulationState-optimalBlendingInv)/(numAges); optBFeature["bInv"] = sum(optimalBlendingInv)/numAges
	features["bAge"] = ((simulationState' * ages)/max(1,sum(simulationState))-optBAge)/optBAge
	features["bInv12"] = sum(simulationState[1:2]-optimalBlendingInv[1:2])/(2); optBFeature["bInv12"] = sum(optimalBlendingInv[1:2])/2
	features["bInv123"] = sum(simulationState[1:3]-optimalBlendingInv[1:3])/(3); optBFeature["bInv123"] = sum(optimalBlendingInv[1:3])/3
	features["bInv1234"] = sum(simulationState[1:4]-optimalBlendingInv[1:4])/(4); optBFeature["bInv1234"] = sum(optimalBlendingInv[1:4])/4
	features["bInv23"] = sum(simulationState[2:3]-optimalBlendingInv[2:3])/(2); optBFeature["bInv23"] = sum(optimalBlendingInv[2:3])/2
	features["bInv234"] = sum(simulationState[2:4]-optimalBlendingInv[2:4])/(3); optBFeature["bInv234"] = sum(optimalBlendingInv[2:4])/3
	features["bInv2345"] = sum(simulationState[2:5]-optimalBlendingInv[2:5])/(4); optBFeature["bInv2345"] = sum(optimalBlendingInv[2:5])/4
	features["bInv34"] = sum(simulationState[3:4]-optimalBlendingInv[3:4])/(2); optBFeature["bInv34"] = sum(optimalBlendingInv[3:4])/2
	features["bInv345"] = sum(simulationState[3:5]-optimalBlendingInv[3:5])/(3); optBFeature["bInv345"] = sum(optimalBlendingInv[3:5])/3
	features["bInv45"] = sum(simulationState[4:5]-optimalBlendingInv[4:5])/(2); optBFeature["bInv45"] = sum(optimalBlendingInv[4:5])/2
	for i in 1:numAges
		features[string("bInv",i)] = simulationState[i]-optimalBlendingInv[i]; optBFeature[string("bInv",i)] = (optimalBlendingInv[i])
	end
	maxFulfillment = Dict()
	for c in wineClasses
		maxFulfillment[c] = maximum(fulfillProduct[a][c] for a in stateBlendingActions[transState])
		features[string("bFFMax",c)] = maxFulfillment[c]/detDemand[c]
	end
	features["bMaxRev"] = maximum(blendingContribution[a] for a in stateBlendingActions[transState])

	#underfulfillment tree
	ffPred = Dict()
	for c in wineClasses
		println("Wine: ", c, "yo")
		uffPrediction = nothing
		foundLeaf = false
		currentNode = "1"
		counter = 0
		while !foundLeaf
			counter += 1
			println(counter)
			if blendingFFTree[c][currentNode]["isLeaf"]
				uffPrediction = blendingFFTree[c][currentNode]["class"]
				println("Rule: ",uffPrediction)
				foundLeaf=true
			else
				currentFeature = blendingFFTree[c][currentNode]["feature"]
				println("Feature: ", currentFeature)
				value = features[currentFeature]
				println("Value of current state: ", value)
				if currentFeature in keys(optBFeature)
					treeValue = blendingFFTree[c][currentNode]["value"] - optBFeature[currentFeature]
				else
					treeValue = blendingFFTree[c][currentNode]["value"]
				end
				println("Value of tree feature: ", treeValue)
				currentNode = value < treeValue ? string(blendingFFTree[c][currentNode]["children"][1]) : string(blendingFFTree[c][currentNode]["children"][2])
			end
		end
		ffPred[c] = Int(round(max(maxFulfillment[c] - uffPrediction,0),digits = 0))
	end
	println("Fulfillment: ", ffPred)

	detPrediction=Dict()
	for c in wineClasses
		println("Wine: ", c, "yo")
		foundLeaf = false
		currentNode = "1"
		if ffPred[c] > 0
			counter = 0
			while !foundLeaf
				counter += 1
				println(counter)
				if blendingDetTree[c][ffPred[c]][currentNode]["isLeaf"]
					detPrediction[c] = blendingDetTree[c][ffPred[c]][currentNode]["class"]
					println("Rule: ",detPrediction[c])
					foundLeaf=true
				else
					currentFeature = blendingDetTree[c][ffPred[c]][currentNode]["feature"]
					println("Feature: ", currentFeature)
					value = features[currentFeature]
					println("Value of current state: ", value)
					if currentFeature in keys(optBFeature)
						treeValue = blendingDetTree[c][ffPred[c]][currentNode]["value"] - optBFeature[currentFeature]
					else
						treeValue = blendingDetTree[c][ffPred[c]][currentNode]["value"]
					end
					println("Value of tree feature: ", treeValue)
					currentNode = value < blendingDetTree[c][ffPred[c]][currentNode]["value"] ? string(blendingDetTree[c][ffPred[c]][currentNode]["children"][1]) : string(blendingDetTree[c][ffPred[c]][currentNode]["children"][2])
				end
			end
		else
			detPrediction[c] = 0
		end
		println("Blending Details: ", detPrediction)

	end
end
#------------------------------------------------#
function isBlend(input)
	if input != []
		return any(input[i]!=input[j] for i in 1:length(input)-1 for j in 2:length(input))
	else
		return false
	end
end
#------------------------------------------------#
function remove!(a, item)
    deleteat!(a, findfirst(x->x==item, a))
end
#-------------------------------------------#
simulationIntegrated()
