using JSON
using DecisionTree
using Distributions
using Dates
using Base.Threads
using Random
#using PyPlot
# using StatsBase
using DataFrames
using Statistics
#using ScikitLearn
#using JLD2
using XLSX
using StatsBase
using Combinatorics
#-------------------------------------------#
#= list of features
Overall inventory level
inventory level >= last target age
inventory level <= first target age
first target age <= inventory level <= last target age
inventory level < first target age
inventory level > last target age
inventory average age
purchase price
number of applicable blending options for older wine
number of applicable blending options for younger wine
number of applicable blending options
maximum demand fulfillment older wine
maximum demand fulfillment younger wine
=#
#-------------------------------------------#
function decisionTree()
	for nP in [2]
		for sF in [25.0/6]
			for pU in [true]
				pSrange = pU ? [20] : ["na"]
				for sprRev in [3]
					for oDP in [0.25]
						for pS in pSrange

							# nP = 2
							# sF = 3
							# sprRev = 3
							# oDP = 0.2
							# pS = 30

							start = now()

							runTrees = true
							findRules = true

							ages = [i for i in 1:6]
							numAges = ages[end]

							scalingFactor = sF

							agesLarge = [i for i in 1:25]
							numAgesLarge = agesLarge[end]

							largeSpread = 10

							ageBuckets = [1,1,1,1,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6]

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
							harvestCosts = nPriceLevels%2 == 0 ? [] : [wineCostMean]
							for i in 1:floor(nPriceLevels/2)
								append!(harvestCosts,wineCostMean - (i-(((nPriceLevels+1)%2)/2)) * costStep)
								append!(harvestCosts,wineCostMean + (i-(((nPriceLevels+1)%2)/2)) * costStep)
							end
							sort!(harvestCosts)
							HtoPC = round(brandContribution[end]/mean(harvestCosts), digits=1)
							println(harvestCosts)

							#maximally allowed spread in blending
							maxSpread = 2

							#indicator for number of harvest price scenarios
							minYield = 0

							#set minimum number of observations per leaf node
							minbucket_size = 0.0001

							possibleMixes = Dict(c=>Dict(i => collect(Set([sort(k) for k in powerset(repeat(copy(ages),outer=[i]),i,i) if sum(k) >= targetAges[c]*i && !(any(any(k1>k2+maxSpread for k2 in k) for k1 in k))])) for i in 1:detDemand[c]) for c in wineClasses)
							#for c in wineClasses for i in 2:detDemand[c] possibleMixes[c][i] = vcat(possibleMixes[c][i],[fill(j,i) for j in targetAges[c]:numAges]) end end
							for c in wineClasses for i in 1:detDemand[c] sort!(possibleMixes[c][i]) end end
							println(possibleMixes)

							#load results from integrated model
							results = Dict()
							currentInstance = ""
							rootComponents=splitpath(@__DIR__)
							for i in 1:length(rootComponents)
								if i==1
									currentInstance *= rootComponents[i][1:2]*"/"
								elseif i==length(rootComponents)
									currentInstance *= rootComponents[i]
								else
									currentInstance *= rootComponents[i]*"/"
								end
							end
							if pU
								simulationInstance = currentInstance * string("/Analysis/Simulation_Aggregated/",nP,"nP_",sF,"sF_","PU","pU_",sprRev,"sprRev_",oDP,"oDP_",pS,"pS_Simulation.json")
								currentInstance *= string("/Results/",nP,"nP_",sF,"sF_","PU","pU_",sprRev,"sprRev_",oDP,"oDP_",pS,"pS_MDP.json")
							else
								simulationInstance = currentInstance * string("/Analysis/Simulation_Aggregated/",nP,"nP_",sF,"sF_","NoPU","pU_",sprRev,"sprRev_",oDP,"oDP_",pS,"pS_Simulation.json")
								currentInstance *= string("/Results/",nP,"nP_",sF,"sF_","NoPU","pU_",sprRev,"sprRev_",oDP,"oDP_",pS,"pS_MDP.json")
							end
							open(currentInstance, "r") do f
								results = JSON.parse(f)
							end

							#load data and classifications
							stateSpaceSize = length(values(results["stateSpace"]))
							purchasingStateSpaceSize = length(values(results["purchasingStateSpace"]))
							stateSpace = Dict(parse(Int32,i)=>results["stateSpace"][i] for i in keys(results["stateSpace"]))
							purchasingStateSpace = Dict(parse(Int32,i)=>results["purchasingStateSpace"][i] for i in keys(results["purchasingStateSpace"]))
							blendingSpace = Dict(parse(Int16,i)=>results["blendingActionSpace"][i] for i in keys(results["blendingActionSpace"]))
							blendingSpaceSize = length(values(blendingSpace))
							fulfillProduct = Dict(a => Dict(w => results["fulfillProduct"][string(a)][string(w)] for w in wineClasses) for a in 1:blendingSpaceSize)
							purchasingSpace = [i for i in 0:nInventoryLevels-1]
							optimalBlending = [results["optimalBlending"][string(s)] for s in 1:stateSpaceSize]
							prodSpecActions = Dict(a => Dict(w => results["prodSpecActions"][string(a)][string(w)] for w in wineClasses) for a in 1:blendingSpaceSize)
							blendingProduct = Dict(a => Dict(w => results["blendingProduct"][string(a)][string(w)] for w in wineClasses) for a in 1:blendingSpaceSize)
							actionPattern = Dict(a => [fulfillProduct[a][w] for w in wineClasses] for a in 1:blendingSpaceSize)
							actionAge = Dict(b => ((blendingSpace[b])' * ages)/max(1,sum(blendingSpace[b])) for b in 1:blendingSpaceSize)
							ffPatternBlending = [actionPattern[optimalBlending[s]] for s in 1:stateSpaceSize]
							ffPatternStateList = Dict(c => Dict(i=>[] for i in 0:detDemand[c]) for c in wineClasses)
							ffPatternDetails = Dict(c => Dict(i=> [] for i in 0:detDemand[c]) for c in wineClasses)
							demandPatterns = sort!(collect(Set(values(actionPattern))), by=x->reverse(x))
							println(demandPatterns)
							blendingDec = Dict(w=>[blendingProduct[optimalBlending[s]][w] for s in 1:stateSpaceSize] for w in wineClasses)
							patternPositions = Dict(i => findfirst(isequal(i),demandPatterns) for i in demandPatterns)
							stateBlendingActions = Dict(parse(Int32,s)=>results["stateBlendingActions"][s] for s in keys(results["stateBlendingActions"]))
							averageAge = Dict(a=>Dict(w=>((prodSpecActions[a][w])' * ages)/max(fulfillProduct[a][w],1) for w in wineClasses) for a in 1:blendingSpaceSize)
							blendingContribution = [sum(fulfillProduct[a][w]*(brandContribution[w] - holdingCosts*(averageAge[a][w]-targetAges[w])) for w in wineClasses) for a in 1:blendingSpaceSize]


							usedStatesBlending = [s for s in 1:stateSpaceSize]
							usedStatesPurchasing = [s for s in 1:purchasingStateSpaceSize]
							ageStructurePurchasing = []
							optimalInvBlending = zeros(numAges)
							optimalInvPurchasing = zeros(numAges-1)

							if isfile(simulationInstance)
								println("using results from simulation instance")
								simulRes = Dict()
								open(simulationInstance, "r") do f
									simulRes = JSON.parse(f)
								end
								stateVisitsBlending = simulRes["visitedBlending"]

								usedStatesBlending = [i for i in 1:stateSpaceSize if stateVisitsBlending[i] >= 1]
								blendingWeights = [stateVisitsBlending[s] for s in usedStatesBlending]
								println("number of used states:")
								println("blending: ",length(usedStatesBlending))
								stateVisitsPurchasing = simulRes["visitedPurchasing"]
								usedStatesPurchasing = [i for i in 1:purchasingStateSpaceSize if stateVisitsPurchasing[i] >= 1]
								purchasingWeights = [stateVisitsPurchasing[s] for s in usedStatesPurchasing]
								println("purchasing: ",length(usedStatesPurchasing))
								ageStructureBlending = [round(simulRes["ageStructureBlending"][i], digits=3) for i in 1:numAges]
								ageStructurePurchasing = [round(simulRes["ageStructurePurchasing"][i], digits=3) for i in 1:numAges]
								optimalAverageAgeB = ((ageStructureBlending)' * ages) / sum(ageStructureBlending)
								optimalAverageAgeP = ((ageStructurePurchasing[1:end-1])' * ages[2:end]) / sum(ageStructurePurchasing[1:end-1])
							else
								println("simulation file doesn't exist")
							end

							#get deliberate underfulfillment decision
							maxBlending = Dict(s => Dict(w => maximum([fulfillProduct[b][w] for b in stateBlendingActions[s]]) for w in wineClasses) for s in 1:stateSpaceSize)
							deliberateUnderfulfillment = Dict(s => Dict(w => maxBlending[s][w] - fulfillProduct[optimalBlending[s]][w] for w in wineClasses) for s in 1:stateSpaceSize)

							#get blending and age details for all
							for s in 1:stateSpaceSize
								for c in wineClasses
									append!(ffPatternStateList[c][fulfillProduct[optimalBlending[s]][c]],s)
									append!(ffPatternDetails[c][fulfillProduct[optimalBlending[s]][c]],[blendingProduct[optimalBlending[s]][c]])
								end
							end

							optimalPurchasing = [results["optimalPurchasing"][string(s)] for s in 1:purchasingStateSpaceSize]
							optimalOrderUpTo = [sum(purchasingStateSpace[s][2:end]) + results["optimalPurchasing"][string(s)] for s in 1:purchasingStateSpaceSize]
							purchaseInvState = [s - purchasingStateSpace[s][1]*Int(purchasingStateSpaceSize/nPriceLevels) for s in 1:purchasingStateSpaceSize]
							averageReward = results["averageReward"]
							results = nothing

							avgDemandAge = sum(detDemand[c] * targetAges[c] for c in wineClasses)/sum(detDemand[c] for c in wineClasses)
							totalDemand = sum(detDemand[c] for c in wineClasses)
							bAge = [round((sum(stateSpace[s][i]*i for i in 1:numAges)/max(1,sum(stateSpace[s]))-optimalAverageAgeB)/optimalAverageAgeB, digits=3) for s in 1:stateSpaceSize]
							bInv = [round(sum(stateSpace[s])/(numAges), digits=3) for s in 1:stateSpaceSize]
							bInv12 = [round(sum(stateSpace[s][i] for i in 1:2)/(2), digits=3) for s in 1:stateSpaceSize]
							bInv123 = [round(sum(stateSpace[s][i] for i in 1:3)/(3), digits=3) for s in 1:stateSpaceSize]
							bInv1234 = [round(sum(stateSpace[s][i] for i in 1:4)/(4), digits=3) for s in 1:stateSpaceSize]
							bInv12345 = [round(sum(stateSpace[s][i] for i in 1:5)/(5), digits=3) for s in 1:stateSpaceSize]
							bInv23 = [round(sum(stateSpace[s][i] for i in 2:3)/(2), digits=3) for s in 1:stateSpaceSize]
							bInv234 = [round(sum(stateSpace[s][i] for i in 2:4)/(3), digits=3) for s in 1:stateSpaceSize]
							bInv2345 = [round(sum(stateSpace[s][i] for i in 2:5)/(4), digits=3) for s in 1:stateSpaceSize]
							bInv23456 = [round(sum(stateSpace[s][i] for i in 2:5)/(5), digits=3) for s in 1:stateSpaceSize]
							bInv34 =  [round(sum(stateSpace[s][i] for i in 3:4)/(2), digits=3) for s in 1:stateSpaceSize]
							bInv345 = [round(sum(stateSpace[s][i] for i in 3:5)/(3), digits=3) for s in 1:stateSpaceSize]
							bInv3456 = [round(sum(stateSpace[s][i] for i in 3:6)/(4), digits=3) for s in 1:stateSpaceSize]
							bInv45 = [round(sum(stateSpace[s][i] for i in 4:5)/(2), digits=3) for s in 1:stateSpaceSize]
							bInv456 = [round(sum(stateSpace[s][i] for i in 4:6)/(3), digits=3) for s in 1:stateSpaceSize]
							bInv56 = [round(sum(stateSpace[s][i] for i in 5:6)/(2), digits=3) for s in 1:stateSpaceSize]
							bInv6 = [stateSpace[s][6] for s in 1:stateSpaceSize]
							bInv5 = [stateSpace[s][5] for s in 1:stateSpaceSize]
							bInv4 = [stateSpace[s][4] for s in 1:stateSpaceSize]
							bInv3 = [stateSpace[s][3] for s in 1:stateSpaceSize]
							bInv2 = [stateSpace[s][2] for s in 1:stateSpaceSize]
							bInv1 = [stateSpace[s][1] for s in 1:stateSpaceSize]
							bMaxFF = Dict(w => [maximum([fulfillProduct[b][w] for b in stateBlendingActions[s]])/detDemand[w] for s in 1:stateSpaceSize] for w in wineClasses)
							bMaxRev = [maximum(blendingContribution[b] for b in stateBlendingActions[s]) for s in 1:stateSpaceSize]

							#purchasing states
							pAge = [round((sum(purchasingStateSpace[s][i]*i for i in 2:numAges)/max(1,sum(purchasingStateSpace[s][2:end])) - optimalAverageAgeP)/optimalAverageAgeP,digits=3) for s in 1:purchasingStateSpaceSize]
							if pU
								pPrice = [harvestCosts[purchasingStateSpace[s][1]+1] for s in 1:purchasingStateSpaceSize]
							end
							pInv = [round(sum(purchasingStateSpace[s][2:end])/(numAges-1), digits=3) for s in 1:purchasingStateSpaceSize]
							pInv23 = [round(sum(purchasingStateSpace[s][i] for i in 2:3)/(2), digits=3) for s in 1:purchasingStateSpaceSize]
							pInv234 = [round(sum(purchasingStateSpace[s][i] for i in 2:4)/(3), digits=3) for s in 1:purchasingStateSpaceSize]
							pInv2345 = [round(sum(purchasingStateSpace[s][i] for i in 2:5)/(4), digits=3) for s in 1:purchasingStateSpaceSize]
							pInv34 =  [round(sum(purchasingStateSpace[s][i] for i in 3:4)/(2), digits=3) for s in 1:purchasingStateSpaceSize]
							pInv345 = [round(sum(purchasingStateSpace[s][i] for i in 3:5)/(3), digits=3) for s in 1:purchasingStateSpaceSize]
							pInv3456 = [round(sum(purchasingStateSpace[s][i] for i in 3:6)/(4), digits=3) for s in 1:purchasingStateSpaceSize]
							pInv45 = [round(sum(purchasingStateSpace[s][i] for i in 4:5)/(2), digits=3) for s in 1:purchasingStateSpaceSize]
							pInv456 = [round(sum(purchasingStateSpace[s][i] for i in 4:6)/(3), digits=3) for s in 1:purchasingStateSpaceSize]
							pInv56 = [round(sum(purchasingStateSpace[s][i] for i in 5:6)/(2), digits=3) for s in 1:purchasingStateSpaceSize]
							pInv6 = [purchasingStateSpace[s][6] for s in 1:purchasingStateSpaceSize]
							pInv5 = [purchasingStateSpace[s][5] for s in 1:purchasingStateSpaceSize]
							pInv4 = [purchasingStateSpace[s][4] for s in 1:purchasingStateSpaceSize]
							pInv3 = [purchasingStateSpace[s][3] for s in 1:purchasingStateSpaceSize]
							pInv2 = [purchasingStateSpace[s][2] for s in 1:purchasingStateSpaceSize]
							pMaxFF = Dict(w => [bMaxFF[w][purchaseInvState[s]] for s in 1:purchasingStateSpaceSize] for w in wineClasses)
							pMaxRev = [bMaxRev[purchaseInvState[s]] for s in 1:purchasingStateSpaceSize]

							#create trees
							println("start creating trees")
							#purchasing tree
							#define labels
							allPurchasingLabels = optimalPurchasing
							purchasingLabels = [optimalPurchasing[s] for s in usedStatesPurchasing]

							#define features
							#purchasing features
							pFeatureList = []
							purchasingFeatures = Array{Float64}(undef,length(usedStatesPurchasing),0)
							pFeatureCounter = 1
							pFeatureNames = Dict()
							#add desired features to the matrix
							purchasingFeatures = hcat(purchasingFeatures,[pInv[s] for s in usedStatesPurchasing]); pFeatureList = vcat(pFeatureList,("pInv")); pFeatureNames[string("x",pFeatureCounter)] = "pInv"; pFeatureCounter += 1
							purchasingFeatures = hcat(purchasingFeatures,[pInv2[s] for s in usedStatesPurchasing]); pFeatureList = vcat(pFeatureList,"pInv2"); pFeatureNames[string("x",pFeatureCounter)] = "pInv2"; pFeatureCounter += 1
							purchasingFeatures = hcat(purchasingFeatures,[pInv3[s] for s in usedStatesPurchasing]); pFeatureList = vcat(pFeatureList,"pInv3"); pFeatureNames[string("x",pFeatureCounter)] = "pInv3"; pFeatureCounter += 1
							purchasingFeatures = hcat(purchasingFeatures,[pInv4[s] for s in usedStatesPurchasing]); pFeatureList = vcat(pFeatureList,"pInv4"); pFeatureNames[string("x",pFeatureCounter)] = "pInv4"; pFeatureCounter += 1
							purchasingFeatures = hcat(purchasingFeatures,[pInv5[s] for s in usedStatesPurchasing]); pFeatureList = vcat(pFeatureList,"pInv5"); pFeatureNames[string("x",pFeatureCounter)] = "pInv5"; pFeatureCounter += 1
							purchasingFeatures = hcat(purchasingFeatures,[pInv6[s] for s in usedStatesPurchasing]); pFeatureList = vcat(pFeatureList,"pInv6"); pFeatureNames[string("x",pFeatureCounter)] = "pInv6"; pFeatureCounter += 1
							purchasingFeatures = hcat(purchasingFeatures,[pInv23[s] for s in usedStatesPurchasing]); pFeatureList = vcat(pFeatureList,"pInv23"); pFeatureNames[string("x",pFeatureCounter)] = "pInv23"; pFeatureCounter += 1
							purchasingFeatures = hcat(purchasingFeatures,[pInv234[s] for s in usedStatesPurchasing]); pFeatureList = vcat(pFeatureList,"pInv234"); pFeatureNames[string("x",pFeatureCounter)] = "pInv234"; pFeatureCounter += 1
							purchasingFeatures = hcat(purchasingFeatures,[pInv2345[s] for s in usedStatesPurchasing]); pFeatureList = vcat(pFeatureList,"pInv2345"); pFeatureNames[string("x",pFeatureCounter)] = "pInv2345"; pFeatureCounter += 1
							purchasingFeatures = hcat(purchasingFeatures,[pInv34[s] for s in usedStatesPurchasing]); pFeatureList = vcat(pFeatureList,"pInv34"); pFeatureNames[string("x",pFeatureCounter)] = "pInv34"; pFeatureCounter += 1
							purchasingFeatures = hcat(purchasingFeatures,[pInv345[s] for s in usedStatesPurchasing]); pFeatureList = vcat(pFeatureList,"pInv345"); pFeatureNames[string("x",pFeatureCounter)] = "pInv345"; pFeatureCounter += 1
							purchasingFeatures = hcat(purchasingFeatures,[pInv3456[s] for s in usedStatesPurchasing]); pFeatureList = vcat(pFeatureList,"pInv3456"); pFeatureNames[string("x",pFeatureCounter)] = "pInv3456"; pFeatureCounter += 1
							purchasingFeatures = hcat(purchasingFeatures,[pInv45[s] for s in usedStatesPurchasing]); pFeatureList = vcat(pFeatureList,"pInv45"); pFeatureNames[string("x",pFeatureCounter)] = "pInv45"; pFeatureCounter += 1
							purchasingFeatures = hcat(purchasingFeatures,[pInv456[s] for s in usedStatesPurchasing]); pFeatureList = vcat(pFeatureList,"pInv456"); pFeatureNames[string("x",pFeatureCounter)] = "pInv456"; pFeatureCounter += 1
							purchasingFeatures = hcat(purchasingFeatures,[pInv56[s] for s in usedStatesPurchasing]); pFeatureList = vcat(pFeatureList,"pInv56"); pFeatureNames[string("x",pFeatureCounter)] = "pInv56"; pFeatureCounter += 1
							purchasingFeatures = hcat(purchasingFeatures,[pAge[s] for s in usedStatesPurchasing]); pFeatureList = vcat(pFeatureList,("pAge")); pFeatureNames[string("x",pFeatureCounter)] = "pAge"; pFeatureCounter += 1
							purchasingFeatures = hcat(purchasingFeatures,[pMaxRev[s] for s in usedStatesPurchasing]); pFeatureList = vcat(pFeatureList,("pMaxRev")); pFeatureNames[string("x",pFeatureCounter)] = "pMaxRev"; pFeatureCounter += 1
							for w in wineClasses
								purchasingFeatures = hcat(purchasingFeatures,[pMaxFF[w][s] for s in usedStatesPurchasing]); pFeatureList = vcat(pFeatureList,(string("pFFMax",w))); pFeatureNames[string("x",pFeatureCounter)] = string("pFFMax",w); pFeatureCounter += 1
							end
							if pU
								purchasingFeatures = hcat(purchasingFeatures,[pPrice[s] for s in usedStatesPurchasing]); pFeatureList = vcat(pFeatureList,("pPrice")); pFeatureNames[string("x",pFeatureCounter)] = "pPrice"; pFeatureCounter += 1
							end

							allPurchasingFeatures = Array{Float64}(undef,purchasingStateSpaceSize,0)
							#add desired features to the matrix
							allPurchasingFeatures = hcat(allPurchasingFeatures,pInv)
							allPurchasingFeatures = hcat(allPurchasingFeatures,pInv2)
							allPurchasingFeatures = hcat(allPurchasingFeatures,pInv3)
							allPurchasingFeatures = hcat(allPurchasingFeatures,pInv4)
							allPurchasingFeatures = hcat(allPurchasingFeatures,pInv5)
							allPurchasingFeatures = hcat(allPurchasingFeatures,pInv6)
							allPurchasingFeatures = hcat(allPurchasingFeatures,pInv23)
							allPurchasingFeatures = hcat(allPurchasingFeatures,pInv234)
							allPurchasingFeatures = hcat(allPurchasingFeatures,pInv2345)
							allPurchasingFeatures = hcat(allPurchasingFeatures,pInv34)
							allPurchasingFeatures = hcat(allPurchasingFeatures,pInv345)
							allPurchasingFeatures = hcat(allPurchasingFeatures,pInv3456)
							allPurchasingFeatures = hcat(allPurchasingFeatures,pInv45)
							allPurchasingFeatures = hcat(allPurchasingFeatures,pInv456)
							allPurchasingFeatures = hcat(allPurchasingFeatures,pInv56)
							allPurchasingFeatures = hcat(allPurchasingFeatures,pAge)
							allPurchasingFeatures = hcat(allPurchasingFeatures,pMaxRev)
							for w in wineClasses
								allPurchasingFeatures = hcat(allPurchasingFeatures,pMaxFF[w])
							end
							if pU
								allPurchasingFeatures = hcat(allPurchasingFeatures,pPrice)
							end

							blendingLabelsFulfillment = Dict(c=>[ffPatternBlending[s][c] for s in usedStatesBlending] for c in wineClasses)
							allBlendingLabelsFulfillment = Dict(c=>[ffPatternBlending[s][c] for s in 1:stateSpaceSize] for c in wineClasses)
							blendingLabelsUnderfulfillment = Dict(c=>[deliberateUnderfulfillment[s][c] for s in usedStatesBlending] for c in wineClasses)
							allBlendingLabelsDetails = ffPatternDetails
							blendingLabelsDetails = Dict(c => Dict(i => [] for i in 1:detDemand[c]) for c in wineClasses)
							ffPatternStateListSimul = Dict(c => Dict(i => [] for i in 1:detDemand[c]) for c in wineClasses)
							for c in wineClasses
								for i in 1:detDemand[c]
									for s in 1:length(ffPatternStateList[c][i])
										if ffPatternStateList[c][i][s] in usedStatesBlending
											append!(ffPatternStateListSimul[c][i],ffPatternStateList[c][i][s])
											append!(blendingLabelsDetails[c][i],[blendingDec[c][ffPatternStateList[c][i][s]]])
										end
									end
								end
							end

							#blending features
							#fulfillment
							bFeatureList = []
							blendingFeatures = Array{Float64}(undef,length(usedStatesBlending),0)
							bFeatureCounter = 1
							bFeatureNames = Dict()
							#add desired features to the matrix
							blendingFeatures = hcat(blendingFeatures,[bInv[s] for s in usedStatesBlending]); bFeatureList = vcat(bFeatureList,("bInv")); bFeatureNames[string("x",bFeatureCounter)] = "bInv"; bFeatureCounter += 1
							blendingFeatures = hcat(blendingFeatures,[bInv1[s] for s in usedStatesBlending]); bFeatureList = vcat(bFeatureList,("bInv1")); bFeatureNames[string("x",bFeatureCounter)] = "bInv1"; bFeatureCounter += 1
							blendingFeatures = hcat(blendingFeatures,[bInv2[s] for s in usedStatesBlending]); bFeatureList = vcat(bFeatureList,("bInv2")); bFeatureNames[string("x",bFeatureCounter)] = "bInv2"; bFeatureCounter += 1
							blendingFeatures = hcat(blendingFeatures,[bInv3[s] for s in usedStatesBlending]); bFeatureList = vcat(bFeatureList,("bInv3")); bFeatureNames[string("x",bFeatureCounter)] = "bInv3"; bFeatureCounter += 1
							blendingFeatures = hcat(blendingFeatures,[bInv4[s] for s in usedStatesBlending]); bFeatureList = vcat(bFeatureList,("bInv4")); bFeatureNames[string("x",bFeatureCounter)] = "bInv4"; bFeatureCounter += 1
							blendingFeatures = hcat(blendingFeatures,[bInv5[s] for s in usedStatesBlending]); bFeatureList = vcat(bFeatureList,("bInv5")); bFeatureNames[string("x",bFeatureCounter)] = "bInv5"; bFeatureCounter += 1
							blendingFeatures = hcat(blendingFeatures,[bInv6[s] for s in usedStatesBlending]); bFeatureList = vcat(bFeatureList,("bInv6")); bFeatureNames[string("x",bFeatureCounter)] = "bInv6"; bFeatureCounter += 1
							blendingFeatures = hcat(blendingFeatures,[bInv12[s] for s in usedStatesBlending]); bFeatureList = vcat(bFeatureList,("bInv12")); bFeatureNames[string("x",bFeatureCounter)] = "bInv12"; bFeatureCounter += 1
							blendingFeatures = hcat(blendingFeatures,[bInv123[s] for s in usedStatesBlending]); bFeatureList = vcat(bFeatureList,("bInv123")); bFeatureNames[string("x",bFeatureCounter)] = "bInv123"; bFeatureCounter += 1
							blendingFeatures = hcat(blendingFeatures,[bInv1234[s] for s in usedStatesBlending]); bFeatureList = vcat(bFeatureList,("bInv1234")); bFeatureNames[string("x",bFeatureCounter)] = "bInv1234"; bFeatureCounter += 1
							blendingFeatures = hcat(blendingFeatures,[bInv12345[s] for s in usedStatesBlending]); bFeatureList = vcat(bFeatureList,("bInv12345")); bFeatureNames[string("x",bFeatureCounter)] = "bInv12345"; bFeatureCounter += 1
							blendingFeatures = hcat(blendingFeatures,[bInv23[s] for s in usedStatesBlending]); bFeatureList = vcat(bFeatureList,("bInv23")); bFeatureNames[string("x",bFeatureCounter)] = "bInv23"; bFeatureCounter += 1
							blendingFeatures = hcat(blendingFeatures,[bInv234[s] for s in usedStatesBlending]); bFeatureList = vcat(bFeatureList,("bInv234")); bFeatureNames[string("x",bFeatureCounter)] = "bInv234"; bFeatureCounter += 1
							blendingFeatures = hcat(blendingFeatures,[bInv2345[s] for s in usedStatesBlending]); bFeatureList = vcat(bFeatureList,("bInv2345")); bFeatureNames[string("x",bFeatureCounter)] = "bInv2345"; bFeatureCounter += 1
							blendingFeatures = hcat(blendingFeatures,[bInv23456[s] for s in usedStatesBlending]); bFeatureList = vcat(bFeatureList,("bInv23456")); bFeatureNames[string("x",bFeatureCounter)] = "bInv23456"; bFeatureCounter += 1
							blendingFeatures = hcat(blendingFeatures,[bInv34[s] for s in usedStatesBlending]); bFeatureList = vcat(bFeatureList,("bInv34")); bFeatureNames[string("x",bFeatureCounter)] = "bInv34"; bFeatureCounter += 1
							blendingFeatures = hcat(blendingFeatures,[bInv345[s] for s in usedStatesBlending]); bFeatureList = vcat(bFeatureList,("bInv345")); bFeatureNames[string("x",bFeatureCounter)] = "bInv345"; bFeatureCounter += 1
							blendingFeatures = hcat(blendingFeatures,[bInv3456[s] for s in usedStatesBlending]); bFeatureList = vcat(bFeatureList,("bInv3456")); bFeatureNames[string("x",bFeatureCounter)] = "bInv3456"; bFeatureCounter += 1
							blendingFeatures = hcat(blendingFeatures,[bInv45[s] for s in usedStatesBlending]); bFeatureList = vcat(bFeatureList,("bInv45")); bFeatureNames[string("x",bFeatureCounter)] = "bInv45"; bFeatureCounter += 1
							blendingFeatures = hcat(blendingFeatures,[bInv456[s] for s in usedStatesBlending]); bFeatureList = vcat(bFeatureList,("bInv456")); bFeatureNames[string("x",bFeatureCounter)] = "bInv456"; bFeatureCounter += 1
							blendingFeatures = hcat(blendingFeatures,[bInv56[s] for s in usedStatesBlending]); bFeatureList = vcat(bFeatureList,("bInv56")); bFeatureNames[string("x",bFeatureCounter)] = "bInv56"; bFeatureCounter += 1
							blendingFeatures = hcat(blendingFeatures,[bAge[s] for s in usedStatesBlending]); bFeatureList = vcat(bFeatureList,("bAge")); bFeatureNames[string("x",bFeatureCounter)] = "bAge"; bFeatureCounter += 1
							blendingFeatures = hcat(blendingFeatures,[bMaxRev[s] for s in usedStatesBlending]); bFeatureList = vcat(bFeatureList,("bMaxRev")); bFeatureNames[string("x",bFeatureCounter)] = "bMaxRev"; bFeatureCounter += 1
							for w in wineClasses
								blendingFeatures = hcat(blendingFeatures,[bMaxFF[w][s] for s in usedStatesBlending]); bFeatureList = vcat(bFeatureList,(string("bFFMax",w))); bFeatureNames[string("x",bFeatureCounter)] = string("bFFMax",w); bFeatureCounter += 1
							end

							#blending features
							allBlendingFeatures = Array{Float64}(undef,stateSpaceSize,0)
							#add desired features to the matrix
							allBlendingFeatures = hcat(allBlendingFeatures,bInv)
							allBlendingFeatures = hcat(allBlendingFeatures,bInv1)
							allBlendingFeatures = hcat(allBlendingFeatures,bInv2)
							allBlendingFeatures = hcat(allBlendingFeatures,bInv3)
							allBlendingFeatures = hcat(allBlendingFeatures,bInv4)
							allBlendingFeatures = hcat(allBlendingFeatures,bInv5)
							allBlendingFeatures = hcat(allBlendingFeatures,bInv6)
							allBlendingFeatures = hcat(allBlendingFeatures,bInv12)
							allBlendingFeatures = hcat(allBlendingFeatures,bInv123)
							allBlendingFeatures = hcat(allBlendingFeatures,bInv1234)
							allBlendingFeatures = hcat(allBlendingFeatures,bInv12345)
							allBlendingFeatures = hcat(allBlendingFeatures,bInv23)
							allBlendingFeatures = hcat(allBlendingFeatures,bInv234)
							allBlendingFeatures = hcat(allBlendingFeatures,bInv2345)
							allBlendingFeatures = hcat(allBlendingFeatures,bInv23456)
							allBlendingFeatures = hcat(allBlendingFeatures,bInv34)
							allBlendingFeatures = hcat(allBlendingFeatures,bInv345)
							allBlendingFeatures = hcat(allBlendingFeatures,bInv3456)
							allBlendingFeatures = hcat(allBlendingFeatures,bInv45)
							allBlendingFeatures = hcat(allBlendingFeatures,bInv456)
							allBlendingFeatures = hcat(allBlendingFeatures,bInv56)
							allBlendingFeatures = hcat(allBlendingFeatures,bAge)
							allBlendingFeatures = hcat(allBlendingFeatures,bMaxRev)
							for w in wineClasses
								allBlendingFeatures = hcat(allBlendingFeatures,bMaxFF[w])
							end


							blendingDetailFeatures = Dict()
							allBlendingDetailFeatures = Dict()
							blendingDetWeights = Dict()
							for c in wineClasses
								blendingDetailFeatures[c] = Dict()
								blendingDetWeights[c] = Dict()
								allBlendingDetailFeatures[c] = Dict()
								for i in 1:detDemand[c]
									blendingDetailFeatures[c][i] = zeros(length(ffPatternStateListSimul[c][i]), length(bFeatureList))
									allBlendingDetailFeatures[c][i] = zeros(length(ffPatternStateList[c][i]), length(bFeatureList))
									blendingDetWeights[c][i] = zeros(length(ffPatternStateListSimul[c][i]))
									for s in 1:length(ffPatternStateListSimul[c][i])
										blendingDetailFeatures[c][i][s,:] = copy(allBlendingFeatures[ffPatternStateListSimul[c][i][s],:])
										blendingDetWeights[c][i][s] = copy(stateVisitsBlending[ffPatternStateListSimul[c][i][s]])
									end
									for s in 1:length(ffPatternStateList[c][i])
										allBlendingDetailFeatures[c][i][s,:] = copy(allBlendingFeatures[ffPatternStateList[c][i][s],:])
									end
								end
							end

							if runTrees
								#use optimal trees by Bertsimas
								println("start Bertsimas trees")
								trees = Dict()
								pLearner = uLearner = dLearner = Dict()
								pPred = nothing
								bPred = Dict()
								dPred = Dict(c=>Dict() for c in wineClasses)
								treeDepthLevels = [6,7,8]
								for tD in treeDepthLevels
									trees[tD] = Dict()
									sampleWeights = Dict(true => [purchasingWeights, blendingWeights, blendingDetWeights], false =>[nothing, nothing, Dict(w=>Dict(i=>nothing for i in 1:detDemand[w]) for w in wineClasses)])
									pLearner[tD] = uLearner[tD] = dLearner[tD] = Dict()
									for sW in [true, false]
										trees[tD][sW] = Dict()
										println("purchasing trees")
										pLearner[tD][sW] = IAI.OptimalTreeClassifier(max_depth=tD, cp=0, minbucket = minbucket_size)
										IAI.fit!(pLearner[tD][sW], purchasingFeatures, purchasingLabels, sample_weight = sampleWeights[sW][1])
										trees[tD][sW]["purchasing"] = getRules(pLearner[tD][sW], pFeatureNames, false)
										if tD == treeDepthLevels[end] && sW == true
											pPred = IAI.predict(pLearner[tD][true], allPurchasingFeatures)
										end
										# p_renamed = IAI.TreePlot(pLearner[tD][sW], feature_renames=pFeatureNames)
										# IAI.show_in_browser(p_renamed)
										# println(IAI.variable_importance(pLearner[tD][sW]))
										# println(IAI.score(pLearner[tD][sW], purchasingFeatures, purchasingLabels))
										println("underfulfillment trees")
										trees[tD][sW]["blendingFF"] = Dict()
										uLearner[tD][sW] = Dict()
										for c in wineClasses
											println(c)
											uLearner[tD][sW][c] = IAI.OptimalTreeClassifier(max_depth=tD, minbucket = minbucket_size, cp=0)
											IAI.fit!(uLearner[tD][sW][c], blendingFeatures, blendingLabelsUnderfulfillment[c], sample_weight = sampleWeights[sW][2])
											if tD == treeDepthLevels[end] && sW == true
												bPred[c] = IAI.predict(uLearner[tD][true][c], allBlendingFeatures)
											end
											trees[tD][sW]["blendingFF"][c] = getRules(uLearner[tD][sW][c], bFeatureNames, false)
											# u_renamed = IAI.TreePlot(uLearner[c], feature_renames=bFeatureNamesFF)
											# IAI.show_in_browser(u_renamed)
											# println(IAI.variable_importance(uLearner[tD][sW][c]))
											# println(IAI.score(uLearner[tD][sW][c], blendingFeaturesFF, blendingLabelsUnderfulfillment[c]))
										end
										dLearner[tD][sW] = Dict(c=> Dict() for c in wineClasses)
										println("blending detail trees")
										for c in wineClasses
											println(c)
											for i in 1:detDemand[c]
												println("demand fulfillment: ", i)
												dLearner[tD][sW][c][i] = IAI.OptimalTreeClassifier(max_depth=tD, minbucket = minbucket_size, cp=0)
												IAI.fit!(dLearner[tD][sW][c][i], blendingDetailFeatures[c][i], blendingLabelsDetails[c][i], sample_weight = sampleWeights[sW][3][c][i])
												trees[tD][sW][string("blendingDet",i,c)] = getRules(dLearner[tD][sW][c][i], bFeatureNames, false)
												if tD == treeDepthLevels[end] && sW == true
													dPred[c][i] = IAI.predict(dLearner[tD][true][c][i], allBlendingFeatures)
												end
												# println(IAI.variable_importance(dLearner[tD][sW][c][i]))
												# if i > 1
												# 	d_renamed = IAI.TreePlot(dLearner[tD][sW][c][i], feature_renames=bFeatureNamesDet[c])
												# 	IAI.show_in_browser(d_renamed)
												# end
												# println(IAI.score(dLearner[tD][sW][c][i], blendingDetailFeatures[c][i], blendingLabelsDetails[c][i]))
											end
										end
									end
								end
								currentInstance = ""
								rootComponents=splitpath(@__DIR__)
								for i in 1:length(rootComponents)
									if i==1
										currentInstance *= rootComponents[i][1:2]*"/"
									elseif i==length(rootComponents)
										currentInstance *= rootComponents[i]
									else
										currentInstance *= rootComponents[i]*"/"
									end
								end
								if pU
									currentInstance *= string("/DecisionTree/PurchasingAmountTrees/",nP,"nP_",sF,"sF_","PU","pU_",sprRev,"sprRev_",oDP,"oDP_",pS,"pS")
								else
									currentInstance *= string("/DecisionTree/PurchasingAmountTrees/",nP,"nP_",sF,"sF_","NoPU","pU_",sprRev,"sprRev_",oDP,"oDP_",pS,"pS")
								end
								trees["pFeatureList"] = pFeatureNames
								trees["bFeatureList"] = bFeatureNames
								open(string(currentInstance, "_Trees.json"), "w") do f
									JSON.print(f, trees, 4)
								end
							end

							#if translate
							results = Dict()

							if findRules
								currentInstance = ""
								rootComponents=splitpath(@__DIR__)
								for i in 1:length(rootComponents)
									if i==1
										currentInstance *= rootComponents[i][1:2]*"/"
									elseif i==length(rootComponents)
										currentInstance *= rootComponents[i]
									else
										currentInstance *= rootComponents[i]*"/"
									end
								end
								if pU
									currentInstance *= string("/DecisionTree/PurchasingAmountTrees/",nP,"nP_",sF,"sF_","PU","pU_",sprRev,"sprRev_",oDP,"oDP_",pS,"pS")
								else
									currentInstance *= string("/DecisionTree/PurchasingAmountTrees/",nP,"nP_",sF,"sF_","NoPU","pU_",sprRev,"sprRev_",oDP,"oDP_",pS,"pS")
								end
								trees = Dict()
								open(string(currentInstance, "_Trees.json"), "r") do f
									trees = JSON.parse(f)
								end
								pTD = 8; pSW = true;
								bfTD = 8; bfSW = true;
								bdTD = Dict(c=>[8 for i in 1:detDemand[c]] for c in wineClasses); bdSW = Dict(c=>[true for i in detDemand[c]] for c in wineClasses);

							end

							purchasingPredictions = [Int(findTreeRule(allPurchasingFeatures[s,:], pFeatureList, trees[string(pTD)][string(pSW)]["purchasing"])) for s in 1:purchasingStateSpaceSize]
							blendingPredictions = Dict(w => [max(0,maxBlending[s][w] - Int(findTreeRule(allBlendingFeatures[s,:], bFeatureList, trees[string(bfTD)][string(bfSW)]["blendingFF"][string(w)]))) for s in 1:stateSpaceSize] for w in wineClasses)#[bPred[s] for s in 1:stateSpaceSize]#

							predictedActions = fill(-1,stateSpaceSize)
							maxAlternatives = 0
							avgAlternatives = 0
							minAlternatives = 1000
							for s in 1:stateSpaceSize
								predictedPattern = vcat([blendingPredictions[w][s] for w in 1:nProducts-1],[maxBlending[s][nProducts]])
								patternPos = findfirst(isequal(predictedPattern),demandPatterns)
								if patternPos == nothing
									println(predictedPattern)
									println(demandPatterns)
								end
								actionAlternatives = []
								while actionAlternatives == []
									predictedPattern = demandPatterns[patternPos]
									for a in stateBlendingActions[s]
										if actionPattern[a] == predictedPattern
											append!(actionAlternatives,a)
										end
									end
									patternPos = patternPos > 1 ? patternPos - 1 : length(demandPatterns)
								end
								blendingPredictionDetails = Dict()
								for c in wineClasses
									if predictedPattern[c] > 0
										blendingPredictionDetails[c] = findTreeRule(allBlendingFeatures[s,:], bFeatureList, trees[string(bdTD[c][Int(predictedPattern[c])])][string(bdSW[c][Int(predictedPattern[c])])][string("blendingDet",Int(predictedPattern[c]),c)])
									else
										blendingPredictionDetails[c] = 0
									end
								end

								#proceed from older to younger wine
								for w in nProducts:-1:1
									predictedDetails = blendingPredictionDetails[w]
									newAlternatives = []
									if predictedDetails != 0
										originalPos = findfirst(isequal(predictedDetails),possibleMixes[w][predictedPattern[w]])
										predictionPos = copy(originalPos)
										counter = 0
										while newAlternatives == []
											counter += 1
											if counter > length(possibleMixes[w][predictedPattern[w]])
												println(stateSpace[s])
												println(predictedDetails)
												error("deadlock")
											end
											for a in actionAlternatives
												if isequal(blendingProduct[a][w],predictedDetails)
													append!(newAlternatives,a)
												end
											end
											if predictionPos <= originalPos + 0.5
												newPos = predictionPos > 1 ? predictionPos - 1 : min(length(possibleMixes[w][predictedPattern[w]]),originalPos + 1)
											else
												newPos = min(length(possibleMixes[w][predictedPattern[w]]),predictionPos+1)
											end
											predictedDetails = possibleMixes[w][predictedPattern[w]][newPos]
											predictionPos = copy(newPos)
										end
										actionAlternatives = copy(newAlternatives)
									end
								end
								finalAlternatives = copy(actionAlternatives)
								maxAlternatives = max(length(finalAlternatives), maxAlternatives)
								minAlternatives = min(length(finalAlternatives), minAlternatives)
								avgAlternatives += length(finalAlternatives)
								if length(finalAlternatives) >= 9
									println(stateSpace[s])
									println(blendingPredictionDetails)
									for a in finalAlternatives
										println(blendingSpace[a])
									end
									error("too many alternatives")
								end
								takenAction = -1
								oldestAlternative = -1
								minDifference = Inf
								for a in finalAlternatives
									if !isempty(ageStructurePurchasing)
										differenceFromOptimal = sum(map(x -> abs(x)^2,(stateSpace[s]-blendingSpace[a])[1:numAges] - ageStructurePurchasing))
										if differenceFromOptimal < minDifference
											minDifference = differenceFromOptimal
											takenAction = a
										end
									else
										if actionAge[a] > oldestAlternative
											oldestAlternative = actionAge[a]
											takenAction = a
										end
									end
								end
								predictedActions[s] = takenAction
							end

							println("avgAlternatives: ", string(avgAlternatives/stateSpaceSize))
							println("maxAlternatives: ", string(maxAlternatives))
							results["purchasingPredictions"] = purchasingPredictions
							results["blendingPredictions"] = predictedActions
							open(string(currentInstance,"_Predictions.json"), "w") do f
								JSON.print(f,results,4)
							end
						end
					end
				end
			end
		end
	end
end
#-------------------------------------------#
function print_tree_to_xlsx(model)
	print("")
end
#-------------------------------------------#
function getHyperRules(learner, features, use_regression)
	tree = Dict()
	for i in 1:IAI.get_num_nodes(learner)
		tree[i] = Dict()
		tree[i]["isLeaf"] = IAI.is_leaf(learner, i)
		if tree[i]["isLeaf"]
			if use_regression
				tree[i]["class"] = Int(round(IAI.get_regression_constant(learner, i), digits=0))
			else
				tree[i]["class"] = IAI.get_classification_label(learner, i)
			end
		else
			tree[i]["hyperplane"] = IAI.is_hyperplane_split(learner, i)
			if IAI.is_hyperplane_split(learner, i)
				tree[i]["children"] = [IAI.get_lower_child(learner, i), IAI.get_upper_child(learner, i)]
				weights = IAI.get_split_weights(learner, i)[1]
				tree[i]["weights"] = Dict()
				for k in keys(weights)
					tree[i]["weights"][features[string(k)]] = weights[k]
				end
				tree[i]["value"] = IAI.get_split_threshold(learner, i)
			else
				tree[i]["children"] = [IAI.get_lower_child(learner, i), IAI.get_upper_child(learner, i)]
				tree[i]["feature"] = features[string(IAI.get_split_feature(learner,i))]
				tree[i]["value"] = IAI.get_split_threshold(learner,i)
			end
		end
	end
	return tree
end
#-------------------------------------------#
function getRules(learner, features, use_regression)
	tree = Dict()
	for i in 1:IAI.get_num_nodes(learner)
		tree[i] = Dict()
		tree[i]["isLeaf"] = IAI.is_leaf(learner, i)
		tree[i]["hyperplane"] = false
		if tree[i]["isLeaf"]
			tree[i]["class"] = IAI.get_classification_label(learner, i)
		else
			tree[i]["children"] = [IAI.get_lower_child(learner, i), IAI.get_upper_child(learner, i)]
			tree[i]["feature"] = features[string(IAI.get_split_feature(learner,i))]
			tree[i]["value"] = IAI.get_split_threshold(learner,i)
		end
	end
	return tree
end
#-------------------------------------------#
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
#-------------------------------------------#
function findTreeRule(features, featureList, tree)
	prediction = nothing
	foundLeaf = false
	currentNode = "1"
	while !foundLeaf
		if tree[currentNode]["isLeaf"]
			prediction = tree[currentNode]["class"]
			foundLeaf=true
		else
			currentFeature = tree[currentNode]["feature"]
			value = features[findfirst(isequal(currentFeature),featureList)]
			treeValue = tree[currentNode]["value"]
			currentNode = value < treeValue ? string(tree[currentNode]["children"][1]) : string(tree[currentNode]["children"][2])
		end
	end
	return prediction
end
#-------------------------------------------#
decisionTree()
