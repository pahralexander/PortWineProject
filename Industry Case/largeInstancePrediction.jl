using JSON
using Distributions
using Dates
using Base.Threads
using Random
#using PyPlot
using StatsBase
using Combinatorics
using XLSX
using DataFrames
#using GLM

#----------------------------------#
function largeInstancePrediction()

		for sF in [25/6.0]
			for nP in [2]
				for pU in [true]
					pSrange = pU ? [20] : ["na"]
					for sprRev in [3]
						for oDP in [0.25]
							for pS in pSrange
								#if !([nP, pU, sprRev, oDP, pS] in [[2,true,3,0.2,15], [2,true,4,0.3,15], [2,true,3,0.3,30]])
								# nP = 2
								# sF = 3
								# sprRev = 3
								# oDP = 0.2
								# pS = 30

								start = now()

								TAp = true
								absLast = true
								useOrderUpTo = false

								agesBase = [i for i in 1:6]
								numAgesBase = agesBase[end]

								scalingFactor = sF

								ages = [i for i in 1:25]
								numAges = ages[end]

								#initialize Wine classes
								nProducts = nP
								wineClasses = [i for i in 1:nProducts]
								targetAgesBase = nProducts > 2 ? [3,4,5] : [3,5]
								targetAges = nProducts > 2 ? [10,15,20] : [10,20]
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
								yieldScenarios = [i for i in 1:nPriceLevels]

								#data representation for saving instance data
								contrExp = sqrt(sprRev)
								contrStart = sprRev == 4 ? 800/3 : 1000/3
								brandContribution = [contrStart*contrExp^(targetAgesBase[i]-targetAgesBase[1]) for i in 1:nProducts]

								#maximally allowed spread in blending
								maxSpread = 2
								largeSpread = 10

								ageBuckets = [1,1,1,1,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6]
								correspondenceSets = Dict(i => [findfirst(isequal(i),ageBuckets),findlast(isequal(i),ageBuckets)] for i in agesBase)
								println(correspondenceSets)

								#initialize data from original large scale case
								supermarketStart = 30/scalingFactor
								supermarketStep = 30/scalingFactor
								maxAgeSupermarket = 3*scalingFactor
								holdingCosts = 30/scalingFactor
								supermarketMax = maxAgeSupermarket*supermarketStep
								supermarketCont = [min(supermarketStart+(i-1)*supermarketStep,supermarketMax) for i in 1:numAges]
								stockValue = [holdingCosts * i for i in 1:numAges]
								supermarketContribution = supermarketCont-stockValue

								#define decay probability in large scale case
								qWeib = 0.875
								βWeib = 0.8
								overallDecay = oDP
								decayProbability = zeros(numAges)
								approxdecayProb = zeros(numAges)
								for k in 1:ages[end]
									approxdecayProb[k] = qWeib^((k-1)^(βWeib))-qWeib^((k)^(βWeib))
								end
								for k in 1:ages[end]
									decayProbability[k] = (approxdecayProb[k]/sum(values(approxdecayProb)))*overallDecay
								end
								println(decayProbability, sum(decayProbability))

								combinedDecay = Dict(s=>zeros(numAges,s+1) for s in 1:nInventoryLevels-1)
								for s in 1:nInventoryLevels-1
									for d in 0:s
										for a in 1:numAges
											combinedDecay[s][a,d+1] = binomial(s,d)*(decayProbability[a]^d)*((1-decayProbability[a])^(s-d))
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

								#indicator for number of harvest price scenarios
								minYield = 0

								if pU
									lambda = exp((2*muH-1)/(2*(sigmaH^2))); qNorm = exp(-1/(sigmaH^2))
								    pmfDivisor=sum((lambda^j)*(qNorm^(j*(j-1)/2)) for j in -200:200)
								    untruncYield = Dict(i=>((lambda^i)*(qNorm^(i*(i-1)/2)))/pmfDivisor for i in minYield:nPriceLevels-1)
								    yieldProbability = [untruncYield[i]/sum(values(untruncYield)) for i in minYield:nPriceLevels-1]
								else
									yieldProbability = [1]
								end
								yieldItems = collect(minYield:nPriceLevels-1)

								maxSalesProd = Dict(c => Dict(i=>0 for i in ages) for c in wineClasses)
								for c in wineClasses
									for i in ages
										if i >= targetAges[c]
											maxSalesProd[c][i] = detDemand[c]
										else
											for k in (detDemand[c]-1):-1:1
												if k*i+(detDemand[c]-k)*min(largeSpread+i,numAges) >= targetAges[c]*detDemand[c]
													maxSalesProd[c][i] = k
												break
												end
											end
										end
									end
								end
								maxSales = Dict()
								for i in ages
									maxSales[i] = sum(maxSalesProd[c][i] for c in wineClasses)
								end

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

								final = now()-start
								println("Time to set up data and probabilities: ",final)

								possibleMixesBase = Dict(c=>Dict(i => sort(collect(Set([sort(k) for k in powerset(repeat(copy(agesBase),inner=[i]),i,i) if sum(k) >= targetAgesBase[c]*i && !(any(any(k1>k2+maxSpread for k2 in k) for k1 in k))])), by=sum) for i in 1:detDemand[c]) for c in wineClasses)
								# println(possibleMixesBase)
								possibleMixes = Dict(c=>Dict(i => sort(collect(Set([sort(k) for k in powerset(repeat(copy(ages),inner=[i]),i,i) if sum(k) >= targetAges[c]*i && !(any(any(k1>k2+largeSpread for k2 in k) for k1 in k))])), by=sum) for i in 1:detDemand[c]) for c in wineClasses)
								# println(possibleMixes)
								targetDeviation = 0
								ageOne = scalingFactor / 2
								spreadLimit = scalingFactor - 1
								baseMixClassification = Dict(c=>Dict(d=>Dict(p=>[] for p in possibleMixes[c][d]) for d in 1:detDemand[c]) for c in wineClasses)
								targetMixes = Dict(c => Dict(d => [p for p in possibleMixesBase[c][d] if mean(p) <= targetAgesBase[c]+targetDeviation] for d in 1:detDemand[c]) for c in wineClasses)
								for c in wineClasses
									for d in 1:detDemand[c]
										for p in possibleMixes[c][d]
											if d==1
												if p[1]-targetAges[c] <= targetDeviation
													baseMixClassification[c][d][p] = [targetAgesBase[c]]
												else
													baseMixClassification[c][d][p] = [max(targetAgesBase[c]+1, ageBuckets[p[1]])]
												end
											else
												spread = p[end]-p[1]
												ageDev = mean(p) - targetAges[c]
												bucketAlloc = [ageBuckets[i] for i in p]
												if ageDev <= targetDeviation && (mean(bucketAlloc) > targetAges[c] || !(bucketAlloc in possibleMixesBase[c][d]))
													minDev = Inf
													for t in targetMixes[c][d]
														if sum(abs.(bucketAlloc-t)) < minDev
														 	baseMixClassification[c][d][p] = t
														 	minDev = sum(abs.(bucketAlloc-t))
														end
													end
												elseif ageDev > targetDeviation
													if mean(bucketAlloc) == targetAgesBase[c]
														if ageDev <= ageOne
															if bucketAlloc[end] > bucketAlloc[1]
																bucketAlloc[1] += 1
															else
																bucketAlloc[end] += 1
															end
														else
															bucketAlloc[1] += 1; bucketAlloc[end] += 1
														end
													else
														counter = 0
														while !(bucketAlloc in possibleMixesBase[c][d])
															bucketAlloc[counter%d + 1] += 1
															counter += 1
														end
													end
													sortedAlloc = sort(bucketAlloc)

													if !(sortedAlloc in possibleMixesBase[c][d])
														allocSum = sum(sortedAlloc)
														allocDev = Inf
														nearestAlloc = possibleMixesBase[c][d][1]
														for m in possibleMixesBase[c][d]
															newDev = sum(m) - allocSum
															if newDev >= 0 && newDev < allocDev
																  nearestAlloc = m
																  allocDev = newDev
															end
														end
														sortedAlloc = nearestAlloc
													end
													baseMixClassification[c][d][p] = sortedAlloc
												else
													baseMixClassification[c][d][p] = sort(bucketAlloc)
												end
											end
										end
									end
								end
								# println(baseMixClassification)

								# load results from simulation and optimization models
								pUString = pU ? "PU" : "NoPU"
								currentInstance = string("",nP,"nP_",sF,"sF_",pUString,"pU_",sprRev,"sprRev_",oDP,"oDP_",pS,"pS")
								println("CURRENT INSTANCE: ",currentInstance)

								optSim = Dict()
								open(string(string(@__DIR__),"/Analysis/Simulation_Aggregated/",currentInstance,"_Simulation.json"), "r") do f
									optSim = JSON.parse(f)
								end
								pricePurchasing = optSim["pricePurchasing"]

								staticOpt = Dict()
								open(string(string(@__DIR__),"/ScaledCase/",currentInstance,"_Opt.json"), "r") do f
									staticOpt = JSON.parse(f)
								end

								pResults = Dict()

								simulationInstance = string("/Analysis/Simulation_Aggregated/",nP,"nP_",sF,"sF_",pUString,"pU_",sprRev,"sprRev_",oDP,"oDP_",pS,"pS_Simulation.json")

								open(string(string(@__DIR__),simulationInstance), "r") do f
									pResults = JSON.parse(f)
								end

								orderUpTo = Int(round(numAges * 2.5, digits = 0))
								orderUpToYoung = Int(round((numAges/2) * 3, digits=0))
								optimalBlendingInv = staticOpt["scaling_approx"]["blendingInv"]
								optimalPurchaseInv = staticOpt["scaling_approx"]["purchasingInv"]
								optimalPurchasingInv = append!(staticOpt["scaling_approx"]["purchasingInv"][2:end], staticOpt["scaling_approx"]["Outdating"])
								optBInvScale = [optimalBlendingInv[correspondenceSets[i][1]] for i in 1:numAgesBase]
								optPInvScale = [optimalPurchaseInv[correspondenceSets[i][1]] for i in 1:numAgesBase]
								println("optimal post-blending inventory: ",optimalPurchasingInv)

								purchasingResults = pResults["priceStockPurchasing"]
								purchasing3Results = pResults["priceStock3Purchasing"]
								rule2DBase = Dict(p=>pResults["2DRule"][string(p)] for p in 0:nPriceLevels-1)

								purchasingStockPrice = Dict(parse(Int,p) => Dict(parse(Int,k)/(numAgesBase-1) => purchasingResults[p][k] == nothing ? -1 : purchasingResults[p][k] for k in keys(purchasingResults[p])) for p in keys(purchasingResults))
								purchasing3StockPrice = Dict(parse(Int,p) => Dict(parse(Int,k)/(numAgesBase-1) => purchasing3Results[p][k] == nothing ? -1 : purchasing3Results[p][k] for k in keys(purchasing3Results[p])) for p in keys(purchasing3Results))
								optBInvBase = pResults["ageStructureBlending"]
								optPInvBase = vcat(pResults["averagePurchasing"],pResults["ageStructurePurchasing"][1:end-1])
								conversionT = sum(optPInvScale[2:end])/sum(optPInvBase[2:end])
								conversion3 = sum(optPInvScale[2:3])/sum(optPInvBase[2:3])
								rule2D = Dict(p=>[Int(round(rule2DBase[p][1]*(sum(optimalPurchaseInv)/sum(optPInvBase)),digits=0)),Int(round(rule2DBase[p][2]*(sum(optimalPurchaseInv[1:correspondenceSets[3][2]])/sum(optPInvBase[1:3])),digits=0))] for p in 0:nPriceLevels-1)
								# rule2D = Dict(p=>[Int(round(rule2DBase[p][1]*sF,digits=0)),Int(round(rule2DBase[p][2]*sF,digits=0))] for p in 0:nPriceLevels-1)
								# rule2D = Dict(p=>[Int(round(rule2DBase[p][1]*sF*(sum(optPInvScale)/sum(optPInvBase)),digits=0)),Int(round(rule2DBase[p][2]*sF*(sum(optPInvScale[1:3])/sum(optPInvBase[1:3])),digits=0))] for p in 0:nPriceLevels-1)

								println("optimal inventory base case: ", optBInvBase, optPInvBase)
								println("expected optimal inventory large case normalized to base case: ", optBInvScale, optPInvScale)
								println("adapted scaling factor: ", string(sum(optPInvScale)/sum(optPInvBase)))
								println("rule2D: ", rule2D)

								optPFeature = Dict("pInv" => sum(optPInvBase[2:end])/5, "pInv2" => optPInvBase[2], "pInv3" => optPInvBase[3], "pInv4" => optPInvBase[4], "pInv5" => optPInvBase[5], "pInv6" => optPInvBase[6], "pInv23" => sum(optPInvBase[2:3])/2,  "pInv234" => sum(optPInvBase[2:4])/3,  "pInv2345" => sum(optPInvBase[2:5])/4, "pInv34" => sum(optPInvBase[3:4])/2, "pInv345" => sum(optPInvBase[3:5])/3, "pInv3456" => sum(optPInvBase[3:6])/4, "pInv45" => sum(optPInvBase[4:5])/2, "pInv456" => sum(optPInvBase[4:6])/3, "pInv56" => sum(optPInvBase[5:6])/2)
								optBFeature = Dict("bInv" => sum(optBInvBase)/6, "bInv1" => optBInvBase[1], "bInv2" => optBInvBase[2], "bInv3" => optBInvBase[3], "bInv4" => optBInvBase[4], "bInv5" => optBInvBase[5], "bInv6" => optBInvBase[6], "bInv12" => sum(optBInvBase[1:2])/2, "bInv123" => sum(optBInvBase[1:3])/3, "bInv1234" => sum(optBInvBase[1:4])/4, "bInv12345" => sum(optBInvBase[1:5])/5, "bInv23" => sum(optBInvBase[2:3])/2, "bInv234" => sum(optBInvBase[2:4])/3, "bInv2345" => sum(optBInvBase[2:5])/4, "bInv23456" => sum(optBInvBase[2:6])/5, "bInv34" => sum(optBInvBase[3:4])/2, "bInv345" => sum(optBInvBase[3:5])/3, "bInv3456" => sum(optBInvBase[3:6])/4, "bInv45" => sum(optBInvBase[4:5])/2, "bInv456" => sum(optBInvBase[4:6])/3, "bInv56" => sum(optBInvBase[5:6])/2)

								#define blending actions
								start=now()
								allActions = defineActionSpace(1,maxSales,numAges,totalDemand,minimumAgeSum,possiblePatterns,targetAges,detDemand,largeSpread,wineClasses)

								blendingSpaceSize = length(allActions[:,1])
								println("Number of actions: ", blendingSpaceSize)
								blendingActionSpace = Dict(a=>allActions[a,1:numAges] for a in 1:blendingSpaceSize)
								stateEncodedBlendings = Dict(a=>(blendingActionSpace[a])' * stateMultipliers +1 for a in 1:blendingSpaceSize)
								demandFulfillment = Dict(a=>allActions[a,numAges+1:numAges+length(wineClasses)] for a in 1:blendingSpaceSize)
								prodSpecActions = Dict(a=>Dict(w=>allActions[a,numAges*w+length(wineClasses)+1:numAges*(w+1)+length(wineClasses)] for w in wineClasses) for a in 1:blendingSpaceSize)

								averageAge = Dict(a => Dict(w => ((prodSpecActions[a][w])' * ages)/max(demandFulfillment[a][w],1) for w in wineClasses) for a in 1:blendingSpaceSize)
								blendingProduct = Dict(a => Dict(w => getBlending(prodSpecActions[a][w]) for w in wineClasses) for a in 1:blendingSpaceSize)
								blendingProductRev = Dict(w => [blendingProduct[a][w] for a in 1:blendingSpaceSize] for w in wineClasses)
								fulfillProduct = Dict(a => Dict(w => demandFulfillment[a][w] for w in wineClasses) for a in 1:blendingSpaceSize)

								actionContribution = Dict(a => sum(fulfillProduct[a][w]*(brandContribution[w] - ((averageAge[a][w]-targetAges[w])*holdingCosts)) for w in wineClasses) for a in 1:blendingSpaceSize); allActions=nothing

								blendingPatterns = Dict(w => collect(Set(copy(blendingProductRev[w]))) for w in wineClasses)
								for w in wineClasses
									remove!(blendingPatterns[w], [])
								end

								blendingClassification = Dict(b=>Dict(c=>[] for c in wineClasses) for b in 1:blendingSpaceSize)
								for b in 1:blendingSpaceSize
									for c in wineClasses
										if length(blendingProduct[b][c]) > 0
											blendingClassification[b][c] = baseMixClassification[c][length(blendingProduct[b][c])][blendingProduct[b][c]]
										end
									end
								end

								#define purchasing actions
								purchasingActionSpace = [i for i in 0:nInventoryLevels-1]

								#load decision trees
								trees = Dict()
								treefile = useOrderUpTo ? string("/DecisionTree/OrderUpToTrees/",currentInstance,"_Trees.json") : string("/DecisionTree/PurchasingAmountTrees/",currentInstance,"_Trees.json")
								open(string(string(@__DIR__),treefile), "r") do f
										trees = JSON.parse(f)
								end
								tDLevels = [6,7,8]
								sWeights = [true, false]
								purchasingTree = Dict(tD => Dict(sW => trees[string(tD)][string(sW)]["purchasing"] for sW in sWeights) for tD in tDLevels)
								blendingFFTree = Dict(tD => Dict(sW => Dict(c=> trees[string(tD)][string(sW)]["blendingFF"][string(c)] for c in wineClasses) for sW in sWeights) for tD in tDLevels)

								pFeatureList = trees["pFeatureList"]
								bFeatureList = trees["bFeatureList"]

								blendingDetTree = Dict(tD => Dict(sW => Dict(c => Dict() for c in wineClasses) for sW in sWeights) for tD in tDLevels)
								for tD in tDLevels
									for sW in sWeights
										for c in wineClasses
											for i in 1:detDemand[c]
												blendingDetTree[tD][sW][c][i] = trees[string(tD)][string(sW)]["blendingDet"*string(i)*string(c)]
											end
										end
									end
								end

								#tree parameter tuning
								allTrees = vcat(["pu"], [string("bf",w) for w in wineClasses if w != wineClasses[end]])

								for w in wineClasses
									for d in 1:detDemand[w]
										allTrees = vcat(allTrees, [string("bd",string(w),string(d))])
									end
								end
								currTD = Dict(t => 8 for t in allTrees); bestTD = Dict(t => 8 for t in allTrees)
								currSW = Dict(); for t in allTrees currSW[t] = false end #t[2] != 'f' end
								bestSW = copy(currSW)
								currConfig = reduce(vcat,[[currTD[t],currSW[t]] for t in allTrees])
								usedConfigs = []

								treeSelectionWeights = [Float64(2^30) for t in allTrees]

								jumpProbability = 0.95

								tDOptions = Dict(t => [6,7,8] for t in allTrees)

								rseeds = [[[StatsBase.sample([i for i in 1:1_000_000], Weights([1/1_000_000 for i in 1:1_000_000])) for i in 1:numAges+1] for it in 1:4000] for k in 1:1]

								tuningStates = Dict(1 => [Int(floor(optimalBlendingInv[i])) for i in 1:numAges], 4 => vcat([Int(floor(optimalBlendingInv[i])) for i in 1:Int(floor(numAges/2))],[Int(ceil(optimalBlendingInv[i])) for i in Int(ceil((numAges+1)/2)):numAges]), 7 =>  vcat([Int(floor(optimalBlendingInv[i]))-1 for i in 1:Int(floor(numAges/2))],[Int(floor(optimalBlendingInv[i])) for i in Int(ceil((numAges+1)/2)):numAges]))
								numTuningRuns = 25

								# if sF > 2
								# 	currentInstanceSmall = string("",nP,"nP_",sF,"sF_",pUString,"pU_",sprRev,"sprRev_",oDP,"oDP_",pS,"pS")
								# 	smallSFFile = string("/Analysis/PurchasingAmount/",currentInstanceSmall,"_ScaleSimulation.json")
								# 	ssfRes = Dict()
								# 	open(string(string(@__DIR__),smallSFFile), "r") do f
								# 			ssfRes = JSON.parse(f)
								# 	end
								# 	currTD = ssfRes["treeDepths"]
								# 	currSW = ssfRes["treeWeights"]
								# 	currConfig = reduce(vcat,[[currTD[t],currSW[t]] for t in allTrees])
								# end

								if TAp == false
									tuningRun = 1
									bestProfits = zeros(10)
									changedTree = changeTD = newFeat = oldFeat = nothing
									while tuningRun <= numTuningRuns
										println("tuning run: ", tuningRun, " test config: ", currConfig)
										sLSeq = Dict(w => [] for w in wineClasses)
										pSeq = []
										seedProfits = []
										append!(usedConfigs, [currConfig]);
										#run experiments with current configuration
										for seedRun in 1:1
											blendingState = tuningStates[seedRun]
											profitSequence = []
											for iteration in 1:1100
												profit = 0
												#apply blending
												appliedBlending, uF = predict_blending(blendingState,blendingFFTree,bFeatureList,blendingDetTree,blendingActionSpace,blendingSpaceSize,demandFulfillment,blendingProduct,blendingClassification,demandPatterns,targetAges,detDemand,numAges,nInventoryLevels,wineClasses,optimalPurchasingInv,optimalBlendingInv,TAp,prodSpecActions,optBInvBase,optBInvScale,scalingFactor,actionContribution, possibleMixesBase, currTD, currSW, absLast, optBFeature, correspondenceSets)
												carryOnState = (blendingState - blendingActionSpace[appliedBlending])[1:numAges-1]
												for w in wineClasses
													append!(sLSeq[w],fulfillProduct[appliedBlending][w])
												end
												profit += actionContribution[appliedBlending]
												#simulate harvest
												Random.seed!(rseeds[seedRun][iteration][1])
												harvestPrice = StatsBase.sample(yieldItems, Weights(yieldProbability))
												#apply purchasing
												purchasingState = vcat([harvestPrice],carryOnState)
												appliedPurchasing = predict_purchasing(Array{Int}(purchasingState),purchasingTree[currTD["pu"]][currSW["pu"]],pFeatureList,blendingActionSpace,blendingSpaceSize,demandFulfillment,demandPatterns,targetAges,detDemand,numAges,nInventoryLevels,wineClasses,harvestCosts,TAp,pricePurchasing,orderUpTo,orderUpToYoung,purchasingStockPrice,purchasing3StockPrice,optimalPurchaseInv,optPInvBase,optPInvScale,actionContribution,scalingFactor,rule2D,useOrderUpTo, pU, optPFeature, correspondenceSets)
												append!(pSeq,appliedPurchasing)
												profit -= (harvestCosts[harvestPrice+1]*appliedPurchasing + holdingCosts*(appliedPurchasing + sum(carryOnState[1:numAges-1])))
												preHarvestState = vcat([appliedPurchasing],carryOnState)
												blendingState = zeros(numAges)
												#simulate decay
												decayProfits = 0; decayedItems = 0;
												for carryOnPos in 1:numAges
													currInv = Int(preHarvestState[carryOnPos])
													if currInv > 0
														Random.seed!(rseeds[seedRun][iteration][carryOnPos + 1])
														decay = StatsBase.sample(decayItems[carryOnPos][currInv], Weights(decayWeights[carryOnPos][currInv]))
														blendingState[carryOnPos] = currInv - decay
														profit += supermarketContribution[carryOnPos]*decay
													end
												end
												if iteration > 100
													append!(profitSequence,profit)
												end
											end
											if tuningRun == 1 && seedRun < 10 && !((seedRun + 1) in keys(tuningStates))
												tuningStates[seedRun+1] = blendingState
											end
											append!(seedProfits,mean(profitSequence))
										end
										sLs = [mean(sLSeq[w]) for w in wineClasses]
										pss = mean(pSeq)
										#update best parameterization
										if mean(seedProfits) > mean(bestProfits)
											bestTD = copy(currTD); bestSW = copy(currSW)
											bestProfits = seedProfits
											println("improvement achieved: new best profit: ", mean(bestProfits), " new best config ", currConfig, " serviceLevels ", sLs, " purchasing ", pss)
										else
											# if tuningRun == 2 || tuningRun == numTuningRuns
											# 	absLast != absLast
											# end
											println("no improvement: ", mean(seedProfits), " serviceLevels ", sLs, " purchasing ", pss)
											currTD = copy(bestTD)
											currSW = copy(bestSW)
										end
										#sample new parameterization
										# if tuningRun == 1 || tuningRun == numTuningRuns - 1
										# 	absLast = !absLast
										# else
										if true
											newConfig = copy(currConfig)
											newTD = copy(currTD)
											newSW = copy(currSW)
											changedTrees = []
											doJump = false
											while newConfig in usedConfigs
												Random.seed!(Random.GLOBAL_RNG)
												doJump = StatsBase.sample([true,false],Weights([jumpProbability, 1-jumpProbability]))
												if doJump
													changedTrees = StatsBase.sample(allTrees, Weights(treeSelectionWeights), Int(floor(length(allTrees)/2)), replace=false)
												else
													changedTrees = [StatsBase.sample(allTrees, Weights(treeSelectionWeights))]
												end
												for ct in changedTrees
													changeTD = StatsBase.sample([true,false], Weights([1,1]))
													if changeTD
														tDopt = setdiff(tDOptions[ct], currTD[ct])
														newTD[ct] = StatsBase.sample(tDopt, Weights([1/length(tDopt) for o in tDopt]))
													else
														newSW[ct] = !currSW[ct]
													end
												end
												newConfig = reduce(vcat,[[newTD[t],newSW[t]] for t in allTrees])
											end
											currTD = copy(newTD); currSW = copy(newSW);
											for ct in changedTrees
												treeSelectionWeights[findfirst(isequal(ct),allTrees)] /= 2
											end
											currConfig = copy(newConfig)
											if doJump
												jumpProbability *= 0.95
											end
										end
										tuningRun += 1
									end
									currTD = bestTD
									currSW = bestSW
									# 	if any(p in [b[1] for b in tabuList[changedTree]] for p in [true, false])
									# 		changeTD = true
									# 		newOptions = [o for o in tDOptions[changedTree] if !(o in [b[1] for b in tabuList[changedTree]]) && o != currTD[changedTree]]
									# 		newFeat = StatsBase.sample(newOptions, Weights([1/length(newOptions) for o in newOptions]))
									# 		oldFeat = currTD[changedTree]
									# 		currTD[changedTree] = newFeat
									# 	else
									# 		Random.seed!(Random.GLOBAL_RNG)
									# 		if length(tabuList[changedTree]) >= 2
									# 			changeTD = false
									# 		else
									# 			changeTD = StatsBase.sample([true,false], Weights([1,1]))
									# 		end
									# 		if changeTD
									# 			newOptions = [o for o in tDOptions[changedTree] if !(o in [b[1] for b in tabuList[changedTree]]) && o != currTD[changedTree]]
									# 			if newOptions == []
									# 				println(changedTree)
									# 				error("no new options")
									# 			end
									# 			Random.seed!(Random.GLOBAL_RNG)
									# 			newFeat = StatsBase.sample(newOptions, Weights([1/length(newOptions) for o in newOptions]))
									# 			oldFeat = currTD[changedTree]
									# 			currTD[changedTree] = newFeat
									# 		else
									# 			newFeat = !currSW[changedTree]
									# 			currSW[changedTree] = newFeat
									# 			oldFeat = !newFeat
									# 		end
									# 	end
									# 	println("changed tree: ",changedTree)
									# 	println("best profits: ",mean(bestProfits))
									# 	println("new feature: ",newFeat)
									# 	println("tabu list: ",tabuList)
									# 	tuningRun += 1
									# 	for t in allTrees
									# 		for k in tabuList[t]
									# 			k[2] -= 1
									# 			if k[2] == 0
									# 				remove!(tabuList[t],k)
									# 			end
									# 		end
									# 	end
									# end
								end

								#choose starting state as average blending state
								blendingState = [Int(floor(optimalBlendingInv[i])) for i in 1:numAges]
								numIterations = 20_000
								warm_up_done = false
								iteration = 0
								blendingFoundIndex = 0
								predictedActionsBlending = Dict()
								predictedActionsPurchasing = Dict()
								predictedStatesBlending = []
								predictedStatesPurchasing = []
								underfulfillmentStates = []
								demandFulfillmentFrequency = Dict(p=>0 for p in demandPatterns)
								blendingStatesAcc = zeros(numAges)
								purchasingStatesAcc = zeros(numAges)
								annualProfits = []
								stateVisitsBlending = Dict()
								stateVisitsPurchasing = Dict()
								blendingFrequency = zeros(blendingSpaceSize)
								purchasingFrequency = Dict(i=>0 for i in 0:nInventoryLevels-1)
								outdatingSequence = []
								start = now()
								println("Start experiments")
								while iteration < numIterations
									if iteration == 1000 && !warm_up_done
										predictedStatesBlending = []
										predictedStatesPurchasing = []
										underfulfillmentStates = []
										demandFulfillmentFrequency = Dict(p=>0 for p in demandPatterns)
										blendingStatesAcc = zeros(numAges)
										purchasingStatesAcc = zeros(numAges)
										annualProfits = []
										stateVisitsBlending = Dict()
										stateVisitsPurchasing = Dict()
										blendingFrequency = zeros(blendingSpaceSize)
										purchasingFrequency = Dict(i=>0 for i in 0:nInventoryLevels-1)
										outdatingSequence = []
										iteration = 0
										warm_up_done = true
										println("warm up done!")
									end

									iteration += 1
									profit = 0

									#apply blending
									if !(blendingState in predictedStatesBlending)
										appliedBlending,underfulfillment = predict_blending(blendingState,blendingFFTree,bFeatureList,blendingDetTree,blendingActionSpace,blendingSpaceSize,demandFulfillment,blendingProduct,blendingClassification,demandPatterns,targetAges,detDemand,numAges,nInventoryLevels,wineClasses,optimalPurchasingInv,optimalBlendingInv,TAp,prodSpecActions,optBInvBase,optBInvScale,scalingFactor,actionContribution, possibleMixesBase, currTD, currSW, absLast, optBFeature, correspondenceSets)
										stateVisitsBlending[blendingState] = 1
										predictedActionsBlending[blendingState] = appliedBlending
										underfulfillment ? append!(underfulfillmentStates,[blendingState]) : nothing
										append!(predictedStatesBlending,[blendingState])
									else
										appliedBlending = predictedActionsBlending[blendingState]
										stateVisitsBlending[blendingState] += 1
										# println("found an already predicted state")
									end
									blendingStatesAcc += blendingState
									blendingFrequency[appliedBlending] += 1
									demandFulfillmentFrequency[demandFulfillment[appliedBlending]] += 1
									outdating = blendingState[end] - blendingActionSpace[appliedBlending][end]
									append!(outdatingSequence,outdating)
									carryOnState = (blendingState - blendingActionSpace[appliedBlending])[1:numAges-1]
									profit += actionContribution[appliedBlending]

									#simulate harvest
									harvestPrice = StatsBase.sample(yieldItems, Weights(yieldProbability))

									purchasingState = vcat([harvestPrice],carryOnState)
									#apply purchasing
									if !(purchasingState in predictedStatesPurchasing)
										appliedPurchasing = predict_purchasing(Array{Int}(purchasingState),purchasingTree[currTD["pu"]][currSW["pu"]],pFeatureList,blendingActionSpace,blendingSpaceSize,demandFulfillment,demandPatterns,targetAges,detDemand,numAges,nInventoryLevels,wineClasses,harvestCosts,TAp,pricePurchasing,orderUpTo,orderUpToYoung,purchasingStockPrice,purchasing3StockPrice,optimalPurchaseInv,optPInvBase,optPInvScale,actionContribution,scalingFactor,rule2D,useOrderUpTo,pU,optPFeature,correspondenceSets)
										predictedActionsPurchasing[purchasingState] = appliedPurchasing
										stateVisitsPurchasing[purchasingState] = 1
										append!(predictedStatesPurchasing,[purchasingState])
									else
										appliedPurchasing = predictedActionsPurchasing[purchasingState]
										stateVisitsPurchasing[purchasingState] += 1
									end
									purchasingStatesAcc += vcat(purchasingState[2:end],[outdating])
									purchasingFrequency[appliedPurchasing] += 1
									profit -= (harvestCosts[harvestPrice+1]*appliedPurchasing + holdingCosts*(appliedPurchasing + sum(carryOnState[1:numAges-1])))
									preHarvestState = vcat([appliedPurchasing],carryOnState)

									blendingState = zeros(numAges)
									#simulate decay
									decayProfits = 0; decayedItems = 0;
									for carryOnPos in 1:numAges
										currInv = preHarvestState[carryOnPos]
										if currInv > 0
											decay = StatsBase.sample(decayItems[carryOnPos][currInv], Weights(decayWeights[carryOnPos][currInv]))
											blendingState[carryOnPos] = currInv - decay
											profit += supermarketContribution[carryOnPos]*decay
										end
									end
									append!(annualProfits,profit)
									if iteration%100 == 0
										println(iteration, " iterations done!")
										if iteration%500 == 0
											println("average reward: ", mean(annualProfits))
											for c in wineClasses
												println("service level ",c,"yo: ",string(sum([sum(demandFulfillmentFrequency[p]*i for p in demandPatterns if p[c]==i) for i in 0:detDemand[c]])/(iteration*detDemand[c])))
											end
											println("average purchasing: ", string(sum(purchasingFrequency[i]*i for i in 0:nInventoryLevels-1)/iteration))
										end
									end
								end
								println(mean(annualProfits))
								avgBlendingInventory = blendingStatesAcc/numIterations
								avgPurchasingInventory = purchasingStatesAcc/numIterations
								mostVisitsBlending = findTopStates(stateVisitsBlending,10)
								# for k in predictedStatesBlending
								# 	if stateVisitsBlending[k] > 1
								# 		println(string(k), " visits: ",string(stateVisitsBlending[k]))
								# 	end
								# end
								mostVisitsPurchasing = findTopStates(stateVisitsPurchasing,10)
								mostBlendingActions = partialsortperm(blendingFrequency,1:10,rev=true)
								averagePurchasing = sum(purchasingFrequency[i]*i for i in 0:nInventoryLevels-1)/numIterations
								averageSelling = [avgBlendingInventory[i]-avgPurchasingInventory[i] for i in 1:numAges]
								averageOutdating = mean(outdatingSequence)
								averageProfits = mean(annualProfits)
								profitsStd = std(annualProfits)
								serviceLevels = Dict()
								for c in wineClasses
									serviceLevels[c] = sum(sum(demandFulfillmentFrequency[p] for p in demandPatterns if p[c]==i)*i for i in 0:detDemand[c]) / detDemand[c]
								end
								results	= Dict("avgPurchasing" => averagePurchasing, "avgSelling" => averageSelling, "avgBlendingInventory" => avgBlendingInventory, "avgPurchasingInventory" => avgPurchasingInventory, "avgSelling" => averageSelling, "avgProfits" => averageProfits, "stdProfits" => profitsStd, "serviceLevels" => serviceLevels, "treeDepths" => currTD, "treeWeights" => currSW)

								if TAp
									results	= Dict("avgPurchasing" => averagePurchasing, "avgSelling" => averageSelling, "avgBlendingInventory" => avgBlendingInventory, "avgPurchasingInventory" => avgPurchasingInventory, "avgSelling" => averageSelling, "avgProfits" => averageProfits, "stdProfits" => profitsStd, "serviceLevels" => serviceLevels)
									jsonPath = useOrderUpTo ? string(string(@__DIR__),"/Analysis/TAp_Results/",currentInstance,"_ScaleSimulation_TAp.json") : string(string(@__DIR__),"/Analysis/TAp_Results/",currentInstance,"_ScaleSimulation_TAp.json")
								else
									results	= Dict("avgPurchasing" => averagePurchasing, "avgSelling" => averageSelling, "avgBlendingInventory" => avgBlendingInventory, "avgPurchasingInventory" => avgPurchasingInventory, "avgSelling" => averageSelling, "avgProfits" => averageProfits, "stdProfits" => profitsStd, "serviceLevels" => serviceLevels, "treeDepths" => currTD, "treeWeights" => currSW)
									jsonPath = useOrderUpTo ? string(string(@__DIR__),"/Analysis/OrderUpTo/",currentInstance,"_ScaleSimulation.json") : string(string(@__DIR__),"/Analysis/PurchasingAmount/",currentInstance,"_ScaleSimulation.json")
								end

								open(jsonPath, "w") do f
									JSON.print(f,results,4)
								end
								#write results to xlsx
								if TAp
									xlsxPath = useOrderUpTo ? string(string(@__DIR__),"/Analysis/OrderUpTo/",currentInstance,"_largeInstance_TAp.xlsx") : string(string(@__DIR__),"/Analysis/PurchasingAmount/",currentInstance,"_largeInstance_TAp.xlsx")
								else
									xlsxPath = useOrderUpTo ? string(string(@__DIR__),"/Analysis/OrderUpTo/",currentInstance,"_largeInstance.xlsx") : string(string(@__DIR__),"/Analysis/PurchasingAmount/",currentInstance,"_largeInstance.xlsx")
								end

								XLSX.openxlsx(xlsxPath, mode="w") do xf
									OVS = xf[1]
									XLSX.rename!(OVS,"largeInstance_Results")

									OVS["A1"] = "Average annual profit"
									OVS["B1"] = averageProfits
									OVS["A2"] = "Profit standard deviation"
									OVS["B2"] = profitsStd
									OVS["E1"] = "number of simulation iterations"
									OVS["F1"] = numIterations

									# OVS["A3"] = "Stock levels purchasing"
									# OVS["B3"] = [i for i in 0:(numAges-1)*(nInventoryLevels-1)]
									# OVS["A4"] = "Iteration count"
									# OVS["B4"] = stockLevel
									# OVS["A5"] = "State count"
									# OVS["B5"] = stockLevelCount
									#
									# OVS["A7"] = "Stock levels blending"
									# OVS["B7"] = [i for i in 0:(numAges)*(nInventoryLevels-1)]
									# OVS["A8"] = "Iteration count"
									# OVS["B8"] = stockLevelBlending
									# OVS["A9"] = "State count"
									# OVS["B9"] = stockLevelBlendingCount

									# OVS["A11"] = "Average age intervals"
									# OVS["B11"] = avgAgeIntervals
									# OVS["A12"] = "Frequency"
									# OVS["B12"] = avgAgeFrequencyPurchasing

									OVS["A14"] = "Age classes"
									OVS["B14"] = ages
									# OVS["A15"] = "Relative use of age in sales: "
									# OVS["B15"] = timesAgeSold/sum(timesAgeSold)
									# OVS["A16"] = "Average selling if available"
									# OVS["B16"] = [mean(ageSellingIfAvailable[a]) for a in 1:numAges]
									OVS["A17"] = "Average inventory level blending"
									OVS["B17"] = avgBlendingInventory
									OVS["A18"] = "Average inventory level purchasing"
									OVS["C18"] = avgPurchasingInventory

									OVS["A25"] = "Demand patterns"
									OVS["B25"] = [string(i) for i in demandPatterns]
									OVS["A26"] = "Frequency"
									OVS["B26"] = [demandFulfillmentFrequency[p] for p in demandPatterns]

									# OVS["A30"] = "Blending paterns 12yo"
									# OVS["B30"] = [string(p) for p in blendingPatterns["young"]]
									# OVS["A31"] = "Frequency"
									# OVS["B31"] = [sum(blendingFrequency[b] for b in blendingPatternActions["young"][p]) for p in blendingPatterns["young"]]
									# OVS["A32"] = "Usage if available"
									# OVS["B32"] = [blendingPatternsFrequency["young"][p]/blendingAvailable["young"][p] for p in blendingPatterns["young"]]
									#
									# OVS["A35"] = "Blending patterns 20yo"
									# OVS["B35"] = [string(p) for p in blendingPatterns["old"]]
									# OVS["A36"] = "Frequency"
									# OVS["B36"] = [sum(blendingFrequency[b] for b in blendingPatternActions["old"][p]) for p in blendingPatterns["old"]]
									# OVS["A37"] = "Usage if available"
									# OVS["B37"] = [blendingPatternsFrequency["old"][p]/blendingAvailable["old"][p] for p in blendingPatterns["old"]]

									OVS["A40"] = "Purchasing decisions"
									OVS["B40"] = [i for i in 0:nInventoryLevels-1]
									OVS["A41"] = "Frequency"
									OVS["B41"] = [purchasingFrequency[i] for i in 0:nInventoryLevels-1]

									OVS["A43"] = "Average outdating"
									OVS["B43"] = averageOutdating

									# OVS["A44"] = "Average no. of decayed items"
									# OVS["B44"] = mean(decayedItemsSequence)

									OVS["A45"] = "visited states blending"
									OVS["B45"] = length(predictedStatesBlending)
									OVS["A46"] = "visited states purchasing"
									OVS["B46"] = length(predictedStatesPurchasing)


									OVS["A47"] = "Demand fulfillment young wine"
									OVS["B47"] = [sum(demandFulfillmentFrequency[p] for p in demandPatterns if p[1]==i) for i in 0:detDemand[1]]
									if nProducts > 2
										OVS["A48"] = "Demand fulfillment medium age wine"
										OVS["B48"] = [sum(demandFulfillmentFrequency[p] for p in demandPatterns if p[2]==i) for i in 0:detDemand[2]]
									end
									OVS["A49"] = "Demand fulfillment old wine"
									OVS["B49"] = [sum(demandFulfillmentFrequency[p] for p in demandPatterns if p[end]==i) for i in 0:detDemand[end]]
									OVS["A50"] = "Deliberate Underfulfillment"
									if underfulfillmentStates != []
										OVS["B50"] = sum(stateVisitsBlending[i] for i in underfulfillmentStates)
									else
										OVS["B50"] = 0
									end

									# OVS["A51"] = "No blending young"
									# OVS["B51"] = sum(blendingFrequency[b] for b in noBlendingActions["young"])
									# # OVS["A52"] = "Percentage if available"
									# # OVS["B52"] = sum(blendingFrequency[b] for b in noBlendingActions["young"])/sum(stateVisitsBlending[s] for s in 1:stateSpaceSize if any(b in noBlendingActions["young"] for b in stateBlendingActions[s]))
									# OVS["A53"] =  "No blending old"
									# OVS["B53"] = sum(blendingFrequency[b] for b in noBlendingActions["old"])
									# OVS["A54"] = "Percentage if available"
									# OVS["B54"] = sum(blendingFrequency[b] for b in noBlendingActions["old"])/sum(stateVisitsBlending[s] for s in 1:stateSpaceSize if any(b in noBlendingActions["old"] for b in stateBlendingActions[s]))

									OVS["A56"]="Most visited blending states"
									OVS["A57"]="Frequency"
									OVS["B56"]=[string(mostVisitsBlending[i][1]) for i in 1:10]
									OVS["B57"]=[mostVisitsBlending[i][2] for i in 1:10]

									OVS["A58"]="Most visited purchasing states"
									OVS["A59"]="Frequency"
									OVS["B58"]=[string(mostVisitsPurchasing[i][1]) for i in 1:10]
									OVS["B59"]=[mostVisitsPurchasing[i][2] for i in 1:10]

									OVS["A60"]="Most frequent blending actions"
									OVS["A61"]="Frequency"
									OVS["B60"]=[string(Array{Int64,1}(blendingActionSpace[b])) for b in mostBlendingActions]
									OVS["B61"]=[blendingFrequency[b] for b in mostBlendingActions]

								end
							#end
							end
						 end
					  end
				   end
			    end
			 end
		  end

#------------------------------------------------#
#recursive function for the configuration of all actions
function defineActionSpace(seed,maxSales,numAges,totalDemand,minimumAgeSum,possiblePatterns,targetAges,detDemand,largeSpread,wineClasses,actionVector=Array{Int8,1}(zeros(numAges*(length(wineClasses)+1)+length(wineClasses))))
	pos = copy(seed)
	if pos >= numAges+1
		check = isReasonable(actionVector,numAges,minimumAgeSum,possiblePatterns,targetAges,detDemand,largeSpread,wineClasses)
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
			return transpose(aV)
		else
			return nothing
		end
	end
	actionSpace = Array{Int8}(undef,0,numAges*(length(wineClasses)+1)+length(wineClasses))
	currentElement = pos
	salesMax = maxSales[currentElement]
	for q in Int8(0):Int8(min(totalDemand,salesMax))
		aV = copy(actionVector)
		aV[pos] = q
		newActions = defineActionSpace(pos+1,maxSales,numAges,totalDemand-q,minimumAgeSum,possiblePatterns,targetAges,detDemand,largeSpread,wineClasses,aV)
		if newActions != nothing
			actionSpace = vcat(actionSpace,newActions)
		end
	end
	return actionSpace
end
#------------------------------------------------#
#checks whether a particular action is valid
function isReasonable(actionVector,numAges,minimumAgeSum,possiblePatterns,targetAges,detDemand,largeSpread,wineClasses)
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
						if p[w] <= 1 || !any(k[i]-k[1] > largeSpread for i in 2:length(k))
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
	occupiedPositions = 0
	blendingPattern = []
	for i in 1:length(action)
		while action[i] > 0
			occupiedPositions += 1
			append!(blendingPattern,i)
			action[i] -= 1
		end
	end
	return blendingPattern
end
#------------------------------------------------#
function predict_blending(startingState,blendingFFTree,bFeatureList,blendingDetTree,blendingActionSpace,blendingSpaceSize,demandFulfillment,blendingProduct,blendingClassification,demandPatterns,targetAges,detDemand,numAges,nInventoryLevels,wineClasses,optimalPurchasingInv,optimalBlendingInv,TAp,prodSpecActions,optBInvBase,optBInvScale,sF,actionContribution, possibleMixesBase,currTD,currSW,absLast,optBFeature,correspondenceSets)
	#get applicable actions
	totalDemand = sum(values(detDemand))
	availableActions = [a for a in 1:blendingSpaceSize if isCompatible(startingState,blendingActionSpace[a],totalDemand)]
	optBAge = ((optimalBlendingInv)' * [i for i in 1:numAges])/sum(optimalBlendingInv)

	ages = [i for i in 1:numAges]
	numAgesBase = Int(numAges/sF)
	maxFulfillment = Dict(c => maximum(demandFulfillment[a][c] for a in availableActions) for c in wineClasses)
	if TAp
		newActions = []
		bestPattern = zeros(length(wineClasses))
		for a in availableActions
			if reverse(demandFulfillment[a]) > reverse(bestPattern)
				newActions = [a]
				bestPattern = demandFulfillment[a]
			elseif isequal(demandFulfillment[a],bestPattern)
				append!(newActions, a)
			end
		end
		availableActions = copy(newActions); newActions = nothing

		for c in wineClasses[end]:-1:1
			if bestPattern[c] > 0
				maxTargetAdh = numAges - targetAges[c]
				remAlternatives = []
				for a in availableActions
					targetAdh = (sum(prodSpecActions[a][c][i]*i for i in ages)/max(1,demandFulfillment[a][c])) - targetAges[c]
					if targetAdh < maxTargetAdh
						maxTargetAdh = targetAdh
						remAlternatives = [a]
					elseif targetAdh == maxTargetAdh
						append!(remAlternatives,a)
					end
				end
				availableActions = copy(remAlternatives); remAlternatives = nothing
			end
		end

		for c in wineClasses[end]:-1:1
			if bestPattern[c] > 1
				bestBlend = blendingProduct[availableActions[1]][c]
				minSpr = numAges
				remAlternatives = []
				for a in availableActions
					blend = blendingProduct[a][c]
					spr = blend[end] - blend[1]
					if spr < minSpr
						minSpr = spr
						bestBlend = blend
						remAlternatives = [a]
					elseif spr == minSpr
						if blend < bestBlend
							bestBlend = blend
							remAlternatives = [a]
						elseif blend == bestBlend
							append!(remAlternatives,a)
						end
					end
				end
				availableActions = copy(remAlternatives); remAlternatives = nothing
			end
		end

		len = length(availableActions)
		if len != 1
			error("$len remaining actions")
		end

		return availableActions[1], false


	# 	oldestAlternative = -1
	# 	minDifference = Inf
	# 	# optimalPurchasingInv = []
	# 	for a in availableActions
	# 		nextState = startingState-blendingActionSpace[a]
	# 		differenceFromOptimal = sum(map(x -> abs(x)^2,nextState - optimalPurchasingInv))
	# 		if differenceFromOptimal < minDifference
	# 			minDifference = differenceFromOptimal
	# 			takenAction = a
	# 		end
	# 	end
	# 	return takenAction,false
	end

	# predict action using tree
	avgDemandAge = sum(detDemand[c] * targetAges[c] for c in wineClasses)/sum(detDemand[c] for c in wineClasses)

	baseInventory = zeros(numAgesBase)
	for i in 1:numAgesBase
		baseInventory[i] = mean(startingState[correspondenceSets[i][1]:correspondenceSets[i][2]])
	end
	if absLast
		baseInventory[numAgesBase] = sum(startingState[correspondenceSets[numAgesBase][1]:correspondenceSets[numAgesBase][2]])
	end

	features = Dict()
	# features["bInv"] = sum(startingState)/(numAges))
	# features["bAge"] = ((startingState' * ages)/max(1,sum(startingState))-optBAge)/optBAge
	# features["bInv12"] = sum(startingState[1:2*sF])/(2*sF)
	# features["bInv123"] = sum(startingState[1:3*sF])/(3*sF)
	# features["bInv1234"] = sum(startingState[1:4*sF])/(4*sF)
	# features["bInv23"] = sum(startingState[sF+1:3*sF])/(2*sF)
	# features["bInv234"] = sum(startingState[sF+1:4*sF])/(3*sF)
	# features["bInv2345"] = sum(startingState[sF+1:5*sF])/(4*sF)
	# features["bInv34"] = sum(startingState[2*sF+1:4*sF])/(2*sF)
	# features["bInv345"] = sum(startingState[2*sF+1:5*sF])/(3*sF)
	# features["bInv45"] = sum(startingState[3*sF+1:5*sF])/(2*sF)
	# for i in 1:numAgesBase
	# 	features[string("bInv",i)] = sum(startingState[(i-1)*sF+1:sF*i])/(sF)
	# end
	# features["bInv"] = sum(baseInventory)/(numAgesBase)
	# features["bAge"] = ((startingState' * ages)/max(1,sum(startingState))-optBAge)/optBAge
	# features["bInv12"] = sum(baseInventory[1:2])/(2)
	# features["bInv123"] = sum(baseInventory[1:3])/(3)
	# features["bInv1234"] = sum(baseInventory[1:4])/(4)
	# features["bInv23"] = sum(baseInventory[2:3])/(2)
	# features["bInv234"] = sum(baseInventory[2:4])/(3)
	# features["bInv2345"] = sum(baseInventory[2:5])/(4)
	# features["bInv34"] = sum(baseInventory[3:4])/(2)
	# features["bInv345"] = sum(baseInventory[3:5])/(3)
	# features["bInv45"] = sum(baseInventory[4:5])/(2)
	# for i in 1:numAgesBase
	# 	features[string("bInv",i)] = baseInventory[i]
	# end
	features["bInv"] = sum(startingState-optimalBlendingInv)/(numAges)
	features["bAge"] = ((startingState' * ages)/max(1,sum(startingState))-optBAge)/optBAge
	features["bInv12"] = sum(startingState[correspondenceSets[1][1]:correspondenceSets[2][end]]-optimalBlendingInv[correspondenceSets[1][1]:correspondenceSets[2][end]])/(correspondenceSets[2][end]-correspondenceSets[1][1]+1)
	features["bInv123"] = sum(startingState[correspondenceSets[1][1]:correspondenceSets[3][end]]-optimalBlendingInv[correspondenceSets[1][1]:correspondenceSets[3][end]])/(correspondenceSets[3][end]-correspondenceSets[1][1]+1)
	features["bInv1234"] = sum(startingState[correspondenceSets[1][1]:correspondenceSets[4][end]]-optimalBlendingInv[correspondenceSets[1][1]:correspondenceSets[4][end]])/(correspondenceSets[4][end]-correspondenceSets[1][1]+1)
	features["bInv12345"] = sum(startingState[correspondenceSets[1][1]:correspondenceSets[5][end]]-optimalBlendingInv[correspondenceSets[1][1]:correspondenceSets[5][end]])/(correspondenceSets[5][end]-correspondenceSets[1][1]+1)
	features["bInv23"] = sum(startingState[correspondenceSets[2][1]:correspondenceSets[3][end]]-optimalBlendingInv[correspondenceSets[2][1]:correspondenceSets[3][end]])/(correspondenceSets[3][end]-correspondenceSets[2][1]+1)
	features["bInv234"] = sum(startingState[correspondenceSets[2][1]:correspondenceSets[4][end]]-optimalBlendingInv[correspondenceSets[2][1]:correspondenceSets[4][end]])/(correspondenceSets[4][end]-correspondenceSets[2][1]+1)
	features["bInv2345"] = sum(startingState[correspondenceSets[2][1]:correspondenceSets[5][end]]-optimalBlendingInv[correspondenceSets[2][1]:correspondenceSets[5][end]])/(correspondenceSets[5][end]-correspondenceSets[2][1]+1)
	features["bInv23456"] = sum(startingState[correspondenceSets[2][1]:correspondenceSets[6][end]]-optimalBlendingInv[correspondenceSets[2][1]:correspondenceSets[6][end]])/(correspondenceSets[6][end]-correspondenceSets[2][1]+1)
	features["bInv34"] = sum(startingState[correspondenceSets[3][1]:correspondenceSets[4][end]]-optimalBlendingInv[correspondenceSets[3][1]:correspondenceSets[4][end]])/(correspondenceSets[4][end]-correspondenceSets[3][1]+1)
	features["bInv345"] = sum(startingState[correspondenceSets[3][1]:correspondenceSets[5][end]]-optimalBlendingInv[correspondenceSets[3][1]:correspondenceSets[5][end]])/(correspondenceSets[5][end]-correspondenceSets[3][1]+1)
	features["bInv3456"] = sum(startingState[correspondenceSets[3][1]:correspondenceSets[6][end]]-optimalBlendingInv[correspondenceSets[3][1]:correspondenceSets[6][end]])/(correspondenceSets[6][end]-correspondenceSets[3][1]+1)
	features["bInv45"] = sum(startingState[correspondenceSets[4][1]:correspondenceSets[5][end]]-optimalBlendingInv[correspondenceSets[4][1]:correspondenceSets[5][end]])/(correspondenceSets[5][end]-correspondenceSets[4][1]+1)
	features["bInv456"] = sum(startingState[correspondenceSets[4][1]:correspondenceSets[6][end]]-optimalBlendingInv[correspondenceSets[4][1]:correspondenceSets[6][end]])/(correspondenceSets[6][end]-correspondenceSets[4][1]+1)
	features["bInv56"] = sum(startingState[correspondenceSets[5][1]:correspondenceSets[6][end]]-optimalBlendingInv[correspondenceSets[5][1]:correspondenceSets[6][end]])/(correspondenceSets[6][end]-correspondenceSets[5][1]+1)
	
	for i in 1:numAgesBase
		features[string("bInv",i)] = sum(startingState[correspondenceSets[i][1]:correspondenceSets[i][end]]-optimalBlendingInv[correspondenceSets[i][1]:correspondenceSets[i][end]])/(correspondenceSets[i][end]-correspondenceSets[i][1]+1)
	end
	if absLast
		features[string("bInv",numAgesBase)] = sum(startingState[correspondenceSets[numAgesBase][1]:correspondenceSets[numAgesBase][end]]-optimalBlendingInv[correspondenceSets[numAgesBase][1]:correspondenceSets[numAgesBase][end]])
	end
	for c in wineClasses
		features[string("bFFMax",c)] = maxFulfillment[c]/detDemand[c]
	end
	features["bMaxRev"] = maximum(actionContribution[a] for a in availableActions)

	ffPred = Dict()
	for c in wineClasses[1:end-1]
		tD = currTD[string("bf",c)]
		sW = currSW[string("bf",c)]
		uffPrediction = nothing
		foundLeaf = false
		currentNode = "1"
		while !foundLeaf
			if blendingFFTree[tD][sW][c][currentNode]["isLeaf"]
				uffPrediction = blendingFFTree[tD][sW][c][currentNode]["class"]
				foundLeaf=true
			else
				currentFeature = blendingFFTree[tD][sW][c][currentNode]["feature"]
				value = features[currentFeature]
				if currentFeature in keys(optBFeature)
					treeValue = blendingFFTree[tD][sW][c][currentNode]["value"] - optBFeature[currentFeature]
				else
					treeValue = blendingFFTree[tD][sW][c][currentNode]["value"]
				end
				currentNode = value < treeValue ? string(blendingFFTree[tD][sW][c][currentNode]["children"][1]) : string(blendingFFTree[tD][sW][c][currentNode]["children"][2])
			end
		end
		ffPred[c] = Int(round(max(maxFulfillment[c] - uffPrediction,0),digits = 0))
	end
	ffPred[wineClasses[end]] = maxFulfillment[wineClasses[end]]
	ffPrediction = [ffPred[c] for c in wineClasses]

	#select action
	actionAlternatives = []
	maxOutdating = nInventoryLevels
	newActions = []
	for a in availableActions
		outdating = (startingState - blendingActionSpace[a])[end]
		if outdating < maxOutdating
			newActions = [a]
			maxOutdating = outdating
		elseif outdating == maxOutdating
			append!(newActions,a)
		end
	end
	availableActions = newActions
	currentPattern = findfirst(isequal(ffPrediction),demandPatterns)
	predictedPattern = demandPatterns[currentPattern]
	if currentPattern == 0
		println(ffPrediction)
	end
	while actionAlternatives == []
		predictedPattern = demandPatterns[currentPattern]
		for a in availableActions
			if demandFulfillment[a] == predictedPattern
				append!(actionAlternatives,a)
			end
		end
		if currentPattern > 1
			currentPattern -= 1
		else
			currentPattern = length(demandPatterns)
		end
	end
	ffPred = Dict(c=>predictedPattern[c] for c in wineClasses)

	#predict details
	detPrediction=Dict()
	for c in wineClasses
		foundLeaf = false
		currentNode = "1"
		if ffPred[c] > 0
			tD = currTD[string("bd",c,ffPred[c])]
			sW = currSW[string("bd",c,ffPred[c])]
			while !foundLeaf
				if blendingDetTree[tD][sW][c][ffPred[c]][currentNode]["isLeaf"]
					detPrediction[c] = blendingDetTree[tD][sW][c][ffPred[c]][currentNode]["class"]
					foundLeaf=true
				else
					currentFeature = blendingDetTree[tD][sW][c][ffPred[c]][currentNode]["feature"]
					value = features[currentFeature]
					if currentFeature in keys(optBFeature)
						treeValue = blendingDetTree[tD][sW][c][ffPred[c]][currentNode]["value"] - optBFeature[currentFeature]
					else
						treeValue = blendingDetTree[tD][sW][c][ffPred[c]][currentNode]["value"]
					end
					currentNode = value < blendingDetTree[tD][sW][c][ffPred[c]][currentNode]["value"] ? string(blendingDetTree[tD][sW][c][ffPred[c]][currentNode]["children"][1]) : string(blendingDetTree[tD][sW][c][ffPred[c]][currentNode]["children"][2])
				end
			end
		else
			detPrediction[c] = 0
		end
	end
	#proceed from old to young
	for c in wineClasses[end]:-1:1
		predictedDetails = detPrediction[c]
		newAlternatives = []
		if predictedDetails != 0
			originalPos = findfirst(isequal(predictedDetails),possibleMixesBase[c][predictedPattern[c]])
			if originalPos == nothing
				println(predictedPattern)
				println(startingState)
				println(predictedDetails)
			end
			predictionPos = copy(originalPos)
			counter = 0
			while newAlternatives == []
				counter += 1
				if counter > length(possibleMixesBase[c][predictedPattern[c]])
					println(startingState)
					println(possibleMixesBase[c][predictedPattern[c]][originalPos])
					println(predictedDetails)
					for a in actionAlternatives
						println(blendingActionSpace[a],blendingClassification[a])
					end
					error("deadlock")
				end
				for a in actionAlternatives
					if isequal(blendingClassification[a][c],predictedDetails)
						append!(newAlternatives,a)
					end
				end
				if predictionPos <= originalPos + 0.5
					newPos = predictionPos > 1 ? predictionPos - 1 : min(length(possibleMixesBase[c][predictedPattern[c]]),originalPos + 1)
				else
					newPos = min(length(possibleMixesBase[c][predictedPattern[c]]),predictionPos+1)
				end
				predictedDetails = possibleMixesBase[c][predictedPattern[c]][newPos]
				predictionPos = copy(newPos)
			end
			actionAlternatives = copy(newAlternatives)
		end
	end

	takenAction = -1
	minDifference = Inf
	# optimalPurchasingInv = []
	for a in actionAlternatives
		nextState = (startingState-blendingActionSpace[a])
		differenceFromOptimal = sum(map(x -> abs(x)^2,nextState - optimalPurchasingInv))
		if differenceFromOptimal < minDifference
			minDifference = differenceFromOptimal
			takenAction = a
		end
	end
	underfulfillment = !any(demandFulfillment[takenAction][c] == maxFulfillment[c] for c in wineClasses)

	if takenAction != -1
		return takenAction,underfulfillment
	else
		error("No action taken")
	end
end
#------------------------------------------------#
function predict_purchasing(startingState,purchasingTree,pFeatureList,blendingActionSpace,blendingSpaceSize,demandFulfillment,demandPatterns,targetAges,detDemand,numAges,nInventoryLevels,wineClasses,harvestCosts,TAp,pricePurchasing,orderUpTo,orderUpToYoung,purchasingStockPrice,purchasing3StockPrice,optimalPurchaseInv,optPInvBase,optPInvScale,actionContribution,sF,rule2D,useOrderUpTo,pU,optPFeature,correspondenceSets)
	# predict action using tree
	totalDemand = sum(values(detDemand))
	numAgesBase = Int(numAges/sF)

	price = harvestCosts[startingState[1]+1]
	priceIndex = startingState[1]

	invState = startingState[correspondenceSets[2][1]:end]
	avgDemandAge = sum(detDemand[c] * targetAges[c] for c in wineClasses)/sum(detDemand[c] for c in wineClasses)
	availableActions = [a for a in 1:blendingSpaceSize if isCompatible(vcat([0],startingState[2:end]),blendingActionSpace[a],totalDemand)]
	maxFulfillment = Dict(c => maximum(demandFulfillment[a][c] for a in availableActions) for c in wineClasses)

	ages = [i for i in correspondenceSets[2][1]:numAges]
	optPAge = ((optimalPurchaseInv[correspondenceSets[2][1]:end])' * ages)/sum(optimalPurchaseInv[correspondenceSets[2][1]:end])

	if TAp
		invLevel = sum(startingState[2:end])
		invLevelNormalized = invLevel * sum(optPInvBase)/sum(optPInvScale)
		invLevelYoung = sum(startingState[2:correspondenceSets[3][2]])
		invCat = round((invLevel/length(invState))*(numAgesBase-1), digits=0)/(numAgesBase-1)
		inv3Cat =  round((invLevelYoung/length(invState))*(numAgesBase-1), digits=0)/(numAgesBase-1)
		#takenAction = StatsBase.sample([i for i in minPurchase:nInventoryLevels-1],Weights(purchasingWeights))
		# pRuleT = Int(round(purchasingStockPrice[startingState[1]][invCat], digits=0))
		# pRule3 = Int(round(purchasing3StockPrice[startingState[1]][inv3Cat], digits=0))
		# if pRuleT < 0 || pRule3 < 0
		# 	if sum(startingState)/numAges > 2.5
		# 		pRule = 0
		# 	else
		# 		pRule = nInventoryLevels-1
		# 	end
		# else
		# 	pRule = max(pRule3,pRuleT)
		# end
		pRule = max(0,rule2D[priceIndex][1]-invLevel,rule2D[priceIndex][2]-invLevelYoung)
		takenAction = min(pRule,nInventoryLevels-1)
		return takenAction
	end

	features = Dict()

	# features["pInv"] = sum(startingState[sF+1:end])/(numAges-sF)
	# for i in 2:numAgesBase
	# 	features[string("pInv",i)] = sum(startingState[sF*(i-1)+1:sF*i])/(sF)
	# end
	# features["pInv23"] = sum(startingState[sF+1:sF*3])/(2*sF)
	# features["pInv234"] = sum(startingState[sF+1:sF*4])/(3*sF)
	# features["pInv34"] = sum(startingState[2*sF+1:sF*4])/(2*sF)
	# features["pInv345"] = sum(startingState[2*sF+1:sF*5])/(3*sF)
	# features["pInv45"] = sum(startingState[3*sF+1:sF*5])/(2*sF)
	features["pInv"] = sum(startingState[2:end]-optimalPurchaseInv[2:end])/(numAges-1)
	features["pInv2"] = sum(startingState[2:correspondenceSets[2][end]]-optimalPurchaseInv[2:correspondenceSets[2][end]])/(correspondenceSets[2][end]-2+1)
	for i in 3:numAgesBase
		features[string("pInv",i)] = sum(startingState[correspondenceSets[i][1]:correspondenceSets[i][end]]-optimalPurchaseInv[correspondenceSets[i][1]:correspondenceSets[i][end]])/(correspondenceSets[i][end]-correspondenceSets[i][1]+1)
	end
	features["pInv23"] = sum(startingState[2:correspondenceSets[3][end]]-optimalPurchaseInv[2:correspondenceSets[3][end]])/(correspondenceSets[3][end]-2+1)
	features["pInv234"] = sum(startingState[2:correspondenceSets[4][end]]-optimalPurchaseInv[2:correspondenceSets[4][end]])/(correspondenceSets[4][end]-2+1)
	features["pInv2345"] = sum(startingState[2:correspondenceSets[5][end]]-optimalPurchaseInv[2:correspondenceSets[5][end]])/(correspondenceSets[5][end]-2+1)
	features["pInv34"] = sum(startingState[correspondenceSets[3][1]:correspondenceSets[4][end]]-optimalPurchaseInv[correspondenceSets[3][1]:correspondenceSets[4][end]])/(correspondenceSets[4][end]-correspondenceSets[3][1]+1)
	features["pInv345"] = sum(startingState[correspondenceSets[3][1]:correspondenceSets[5][end]]-optimalPurchaseInv[correspondenceSets[3][1]:correspondenceSets[5][end]])/(correspondenceSets[5][end]-correspondenceSets[3][1]+1)
	features["pInv3456"] = sum(startingState[correspondenceSets[3][1]:correspondenceSets[6][end]]-optimalPurchaseInv[correspondenceSets[3][1]:correspondenceSets[6][end]])/(correspondenceSets[6][end]-correspondenceSets[3][1]+1)
	features["pInv45"] = sum(startingState[correspondenceSets[4][1]:correspondenceSets[5][end]]-optimalPurchaseInv[correspondenceSets[4][1]:correspondenceSets[5][end]])/(correspondenceSets[5][end]-correspondenceSets[4][1]+1)
	features["pInv456"] = sum(startingState[correspondenceSets[4][1]:correspondenceSets[6][end]]-optimalPurchaseInv[correspondenceSets[4][1]:correspondenceSets[6][end]])/(correspondenceSets[6][end]-correspondenceSets[4][1]+1)
	features["pInv56"] = sum(startingState[correspondenceSets[5][1]:correspondenceSets[6][end]]-optimalPurchaseInv[correspondenceSets[5][1]:correspondenceSets[6][end]])/(correspondenceSets[6][end]-correspondenceSets[5][1]+1)
	
	for c in wineClasses
		features[string("pFFMax",c)] = maxFulfillment[c]/detDemand[c]
	end
	features["pMaxRev"] = maximum(actionContribution[a] for a in availableActions)
	features["pAge"] = ((invState' * ages)/max(1,sum(invState))-optPAge)/optPAge
	features["pPrice"] = price

	#predict purchasing
	pPrediction = nothing
	foundLeaf = false
	currentNode = "1"
	# while !foundLeaf
	# 	if purchasingTree[currentNode]["isLeaf"]
	# 		pPrediction = Int(round(purchasingTree[currentNode]["class"],digits=0))
	# 		foundLeaf=true
	# 	else
	# 		value = features[purchasingTree[currentNode]["feature"]]
	# 		currentNode = value < purchasingTree[currentNode]["value"] ? string(purchasingTree[currentNode]["children"][1]) : string(purchasingTree[currentNode]["children"][2])
	# 	end
	# end
	while !foundLeaf
		if purchasingTree[currentNode]["isLeaf"]
			pPrediction = Int(round(purchasingTree[currentNode]["class"],digits=0))
			foundLeaf=true
		else
			currentFeature = purchasingTree[currentNode]["feature"]
			value = features[currentFeature]
			if currentFeature in keys(optPFeature)#["pInv", "pInv2", "pInv3", "pInv4", "pInv5", "pInv23", "pInv34", "pInv45", "pInv234", "pInv345"]
				treeValue = purchasingTree[currentNode]["value"] - optPFeature[currentFeature]
			else
				treeValue = purchasingTree[currentNode]["value"]
			end
			currentNode = value < treeValue ? string(purchasingTree[currentNode]["children"][1]) : string(purchasingTree[currentNode]["children"][2])
		end
	end

	if pPrediction != nothing
		if useOrderUpTo
			if pU
				return max(0,min(pPrediction*sF - sum(startingState[2:end]),nInventoryLevels-1))
			else
				return max(0,min(Int(round(pPrediction*(sum(optimalPurchaseInv)/sum(optPInvBase)),digits=0)) - sum(startingState[2:end]),nInventoryLevels-1))
			end
		else
			# minAvgDev = Inf
			# newAvgDev = abs((0*sF + (startingState[2:sF])' * [sF-i for i in 1:sF-1])/((sF*(sF+1))/2) - pPrediction)
			# # newAvgDev = abs((0 + sum(startingState[2:sF]))/sF - pPrediction)
			# bestAverage = 0
			# currP = 0
			# while newAvgDev < minAvgDev && currP < nInventoryLevels
			# 	minAvgDev = copy(newAvgDev)
			# 	bestAverage = copy(currP)
			# 	currP += 1
			# 	newAvgDev = abs((currP*sF + (startingState[2:sF])' * [sF-i for i in 1:sF-1])/((sF*(sF+1))/2) - pPrediction)
			# 	# newAvgDev = abs((currP + sum(startingState[2:sF]))/sF - pPrediction)
			# end
			# return bestAverage
			return pPrediction
		end
	else
		error("no purchasing action found")
	end
end
#------------------------------------------------#
function remove!(a, item)
    deleteat!(a, findfirst(x->x==item, a))
end
#------------------------------------------------#
function getState(stateNumber, stateMultipliers, numAges)
	state = []
	invIndicator = copy(stateNumber)
	for i in 1:numAges
		append!(state,Int64(floor(invIndicator/stateMultipliers[i])))
		invIndicator = invIndicator % stateMultipliers[i]
	end
	return state
end
#------------------------------------------------#
function findTopStates(stateVisits,N)
	keyList = []
	valueList = []
	minValue = Inf
	minPos = 0
	for k in keys(stateVisits)
		value = stateVisits[k]
		if length(keyList) < N
			append!(keyList,[k])
			append!(valueList,value)
			if value < minValue
				minPos = length(keyList)
				minValue = value
			end
		else
			if value > minValue
				minPos = findfirst(isequal(minValue),valueList)
				replace!(keyList,keyList[minPos]=>k)
				replace!(valueList,valueList[minPos]=>value,count=1)
				minPos = argmin(valueList)
				minValue = minimum(valueList)
			end
		end
	end
	resultDict = Dict(keyList[i] => valueList[i] for i in 1:N)
	return sort(collect(resultDict), by=x->x[2])
end

#------------------------------------------------#
largeInstancePrediction()
