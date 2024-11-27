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
using GLM
using JuMP
using Gurobi

#----------------------------------#
function optimizeAverage()
	for nP in [2]
		for sF in [3,5]
			for pU in [true]
				pSrange = pU ? [15,20,25,30] : ["na"]
				for sprRev in [2,3,4,5]
					for oDP in [0.25]
						for pS in pSrange

							# nP = 2
							# sF = 3
							# sprRev = 3
							# oDP = 0.2
							# pS = 30
							fixPurchasing = false
							start = now()

							ages = [i for i in 1:5]
							numAges = ages[end]

							scalingFactor = sF

							agesLarge = [i for i in 1:numAges*scalingFactor]
							numAgesLarge = agesLarge[end]

							#maximally allowed spread in blending
							maxSpread = 2
							largeSpread = 3*scalingFactor-1

							ageBuckets = []
							for i in 1:numAges
								ageBuckets = vcat(ageBuckets, fill(i,scalingFactor))
							end
							correspondenceSets = Dict(i=>[k for k in agesLarge if ageBuckets[k]==i] for i in 1:numAges)

							#approximation case
							BAp = false
							SAp = false; protectionLevel = 6
							PrAp = false
							PuAp = false
							TAp = false


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
							qWeib = scalingFactor == 3 ? 0.8415 : 0.8917
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
							targetAges = nProducts > 2 ? [2,3,4] : [2,4]
							targetAgesLarge = [Int(mean([(targetAges[c]-1)*scalingFactor + 1, targetAges[c] * scalingFactor])) for c in wineClasses]
							detDemand = nP > 2 ? [2 for i in 1:nProducts] : [4 for i in 1:nProducts]
							totalDemand = sum(detDemand)

							#initialize Inventory levels
							nInventoryLevels = Int(ceil(totalDemand * 1.5)) + 1
							if pU
								nPriceLevels = 9
							else
								nPriceLevels = 1
							end
							stateMultipliers = [nInventoryLevels^(numAges-i) for i in 1:numAges]
							yieldScenarios = [i for i in 1:nPriceLevels]

							#data representation for saving instance data
							contrExp = sqrt(sprRev)
							start_spr = Dict(2=>4000.0/9, 3=>1000.0/3, 4=>800.0/3, 5=>2000.0/9)
							contrStart = start_spr[sprRev]
							brandContribution = [contrStart*contrExp^(targetAges[i]-targetAges[1]) for i in 1:nProducts]
							holdingCosts = holdingCostsLarge * scalingFactor

							HtoM = round(brandContribution[end]/brandContribution[1], digits=1)
							HtoS = round(brandContribution[end]/supermarketContribution[end], digits=1)
							HtoHC = round(brandContribution[end]/holdingCosts, digits=1)

							sigma_level = Dict(15=>1.125395, 20=>1.52036, 25=>2.04058, 30=>3.10915)
							muH = (nPriceLevels-1)/2
							sigmaH = sigma_level[pS]

							#wine costs (=purchasing costs) of each respective age class
							wineCostMean = 150
							costStep = 20
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
							    yieldProbability = Dict(i=>untruncYield[i]/sum(values(untruncYield)) for i in minYield:nPriceLevels-1)
							else
								yieldProbability = Dict(0=>1)
							end


						    targetBlendRange = Dict(c=>[i for i in targetAges[c]-Int(floor(maxSpread/2)):targetAges[c]+Int(floor(maxSpread/2))] for c in wineClasses)
							targetBlendRangeLarge = Dict(c=>[i for i in targetAgesLarge[c]-Int(floor(largeSpread/2)):targetAgesLarge[c]+Int(floor(largeSpread/2))] for c in wineClasses)

							maxSales = Dict(c => Dict(i=>0 for i in agesLarge) for c in wineClasses)
							for c in wineClasses
								for i in agesLarge
									if i >= targetAgesLarge[c]
										maxSales[c][i] = detDemand[c]
									else
										for k in (detDemand[c]-1):-1:1
											if k*i+(detDemand[c]-k)*min(largeSpread+i,numAgesLarge) >= targetAgesLarge[c]*detDemand[c]
												maxSales[c][i] = k
											break
											end
										end
									end
								end
							end
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

							minAgeBlendingUB = Dict(c=>max(1,targetAges[c]-maxSpread+1) for c in wineClasses)

							minAgeBlendingLargeScale = Dict(c=>max(1,targetAgesLarge[c]-Int(floor(largeSpread/2))) for c in wineClasses)
							minAgeBlendingLargeUB = Dict(c=>max(1,targetAgesLarge[c]-largeSpread+1) for c in wineClasses)

						    #load optimal sales per age class from
							if pU
								currentInstance = string(nP,"nP_",sF,"sF_","PU","pU_",sprRev,"sprRev_",oDP,"oDP_",pS,"pS")
							else
								currentInstance = string(nP,"nP_",sF,"sF_","NoPU","pU_",sprRev,"sprRev_",oDP,"oDP_",pS,"pS")
							end
							avgRes = Dict()

							print("CURRENT INSTANCE: ", currentInstance)

						    open(string(string(@__DIR__),"/Analysis/Simulation_Aggregated/",currentInstance,"_Simulation.json"), "r") do f
						    	avgRes = JSON.parse(f)
						    end

							if fixPurchasing
								scaleRes = Dict()
								open(string(string(@__DIR__),"/Analysis/",currentInstance,"_ScaleSimulation.json"), "r") do f
								   scaleRes = JSON.parse(f)
							   end
							   	averageSalesAge = sum((scaleRes["avgBlendingInventory"][i]-scaleRes["avgPurchasingInventory"][i])*i for i in 1:numAgesLarge)/ sum((scaleRes["avgBlendingInventory"][i]-scaleRes["avgPurchasingInventory"][i]) for i in 1:numAgesLarge)
								sumSales = sum((scaleRes["avgBlendingInventory"][i]-scaleRes["avgPurchasingInventory"][i]) for i in 1:numAgesLarge)
								#println(averageSalesAge)
								averageSalesTarget = sum((scaleRes["serviceLevels"][string(w)]/sum((scaleRes["serviceLevels"][string(w)]) for w in wineClasses))*targetAgesLarge[w] for w in wineClasses)
								#println(averageSalesTarget)
								purchasingCostsContribution = sum((scaleRes["serviceLevels"][string(w)]/10000)*brandContribution[w] for w in wineClasses) - holdingCostsLarge * (averageSalesAge-averageSalesTarget) * sumSales - scaleRes["avgProfits"] - holdingCostsLarge*(scaleRes["avgPurchasing"] + sum(scaleRes["avgPurchasingInventory"])) - sum(supermarketContributionLarge[i]*scaleRes["avgPurchasingInventory"][i-1]*decayProbabilityLarge[i] for i in 2:numAges) - supermarketContributionLarge[1]*decayProbability[1]*scaleRes["avgPurchasing"]

							else
								purchasingCostsContribution = 0
							end
							avgBlendingInv = avgRes["ageStructureBlending"]
							avgSelling = avgRes["averageSelling"]
						    avgPurchasing = avgRes["averagePurchasing"]
							avgServiceLevels = avgRes["serviceLevels"]
							avgAgeBase = ((avgBlendingInv)' * ages)/sum(avgBlendingInv)
							println(avgAgeBase)
							
							avgAgeBucket = round(Int, avgAgeBase)
							minAvgAgeScaled = avgAgeBucket-0.5+(avgAgeBase-minimum(correspondenceSets[avgAgeBucket])+0.5)/length(correspondenceSets[avgAgeBucket])
							println(minAvgAgeScaled)
							
							avgPricePurchasing = Dict(p => (avgRes["pricePurchasing"][string(p)])' * [i for i in 0:nInventoryLevels-1] for p in 0:nPriceLevels-1)
						    intercept = Dict(i=>0.0 for i in 1:numAges-1)
						    slope = Dict(i=>0.0 for i in 1:numAges-1)
						    for i in 1:numAges-1
						        slope[i] = avgSelling[i+1]-avgSelling[i]
						        intercept[i] = avgSelling[i]-slope[i]*i
						    end
						    println(slope, intercept)
						    scaledAges = [ageBuckets[i]-0.5+(i-minimum(correspondenceSets[ageBuckets[i]])+0.5)/length(correspondenceSets[ageBuckets[i]]) for i in 1:numAgesLarge]
						    scaledSelling = [max(0,(intercept[min(max(Int(floor(i)),1),numAges-1)] + slope[min(max(Int(floor(i)),1),numAges-1)] * i)/(numAgesLarge/numAges)) for i in scaledAges]
						    println(scaledAges, scaledSelling)

						    #specify allowed deviation from sales amount
						    allowedFluctuation = 0.1

							ub_base = run_optimization(wineClasses, nInventoryLevels, ages, maxSpread, yieldScenarios, yieldProbability, minAgeBlendingUB, pU, brandContribution, holdingCosts, harvestCosts, targetAges, supermarketContribution, decayProbability, detDemand, targetBlendRange, allowedFluctuation, scaledSelling, maxSales, purchasingCostsContribution, fixPurchasing, avgPricePurchasing, avgServiceLevels, avgBlendingInv, avgSelling, targetAges, minAvgAgeScaled,  correspondenceSets, false)
							ub_large = run_optimization(wineClasses, nInventoryLevels, agesLarge, largeSpread, yieldScenarios, yieldProbability, minAgeBlendingLargeUB, pU, brandContribution, holdingCostsLarge, harvestCosts, targetAgesLarge, supermarketContributionLarge, decayProbabilityLarge, detDemand, targetBlendRangeLarge, allowedFluctuation, scaledSelling, maxSales, purchasingCostsContribution, fixPurchasing, avgPricePurchasing, avgServiceLevels, avgBlendingInv, avgSelling, targetAges, minAvgAgeScaled,  correspondenceSets, false)
							scale_inv_large = run_optimization(wineClasses, nInventoryLevels, agesLarge, largeSpread, yieldScenarios, yieldProbability, minAgeBlendingLargeScale, pU, brandContribution, holdingCostsLarge, harvestCosts, targetAgesLarge, supermarketContributionLarge, decayProbabilityLarge, detDemand, targetBlendRangeLarge, allowedFluctuation, scaledSelling, maxSales, purchasingCostsContribution, fixPurchasing, avgPricePurchasing, avgServiceLevels, avgBlendingInv, avgSelling, targetAges, minAvgAgeScaled, correspondenceSets,  true)

							if fixPurchasing
								open(string(string(@__DIR__),"/ScaledCase/",currentInstance,"_Opt.json"), "r") do f
									results = JSON.parse(f)
								end
								results["ub_fix_p"] = ub_large
								open(string(string(@__DIR__),"/ScaledCase/",currentInstance,"_Opt.json"), "w") do f
									JSON.print(f,results,4)
								end
							else
								results = Dict()
								# open(string(string(@__DIR__),"/ScaledCase/",currentInstance,"_Opt.json"), "r") do f
								# 	results = JSON.parse(f)
								# end
								results["scaling_approx"] = scale_inv_large
								results["ub_base"] =  ub_base
								results["ub_large"] =  ub_large
								open(string(string(@__DIR__),"/ScaledCase/",currentInstance,"_Opt.json"), "w") do f
									JSON.print(f,results,4)
								end
							end
						end
					end
				end
			end
		end
	end
end
#----------------------------------#
function run_optimization(wineClasses, nInventoryLevels, ages, maxSpread, yieldScenarios, yieldProbability, minAgeBlending, pU, brandContribution, holdingCosts, harvestCosts, targetAges, supermarketContribution, decayProbability, detDemand, targetBlendRange, allowedFluctuation, scaledSelling, maxSales, purchasingCostContribution, fixPurchasing, avgPricePurchasing, avgServiceLevels, avgBlendingInv, avgSelling, targetAgesBase, minAvgAgeScaled, correspondenceSets, scalePostBlending)

	numAges = ages[end]

	blendingPatterns = Dict(c=>[] for c in wineClasses)
	if scalePostBlending
		blendingCandidates = Dict(c=> [k for k in powerset(ages[minAgeBlending[c]:end],1,detDemand[c])] for c in wineClasses)
	else
		blendingCandidates = Dict(c=> [k for k in powerset(ages[minAgeBlending[c]:end],1,detDemand[c])] for c in wineClasses)
	end
	for c in wineClasses
		for k in blendingCandidates[c]
			if length(k) == 1
				if k[1] >= targetAges[c]
					append!(blendingPatterns[c],[k])
				end
			else
				if maximum(k) > targetAges[c] && maximum(k) - minimum(k) <= maxSpread
						append!(blendingPatterns[c],[k])
				end
			end
		end
	end

	#gurobi model
	model = Model(Gurobi.Optimizer)

	#variable definitions
	@variable(model, usage[c=wineClasses, k=blendingPatterns[c], [i for i in k]], lower_bound=0)
	@variable(model, sales[ages], lower_bound=0)
	@variable(model, purchasing[yieldScenarios], lower_bound=0, upper_bound=nInventoryLevels-1)
	@variable(model, blendingInv[ages], lower_bound=0)
	@variable(model, purchasingInv[ages], lower_bound=0)
	@variable(model, salesDiv[c=wineClasses, minAgeBlending[c]:numAges], lower_bound=0)
	@variable(model, outdating, lower_bound=0)
	@variable(model, factors[1:5])
	@variable(model, constant)

	#objective
	if pU
		obj = @expression(model, sum(sum(salesDiv[c,i] for i in minAgeBlending[c]:numAges)*brandContribution[c] for c in wineClasses) - holdingCosts*sum((sum(salesDiv[c,i]*i for i in minAgeBlending[c]:numAges) - sum(salesDiv[c,i] for i in minAgeBlending[c]:numAges)*targetAges[c]) for c in wineClasses) + sum(purchasingInv[i]*decayProbability[i]*supermarketContribution[i] for i in ages) + outdating*supermarketContribution[end] - holdingCosts*sum(purchasingInv[i] for i in ages) - sum(purchasing[h]*harvestCosts[h]*yieldProbability[h-1] for h in yieldScenarios))
	else
		obj = @expression(model, sum(sum(salesDiv[c,i] for i in minAgeBlending[c]:numAges)*brandContribution[c] for c in wineClasses) - holdingCosts*sum((sum(salesDiv[c,i]*i for i in minAgeBlending[c]:numAges) - sum(salesDiv[c,i] for i in minAgeBlending[c]:numAges)*targetAges[c]) for c in wineClasses) + sum(purchasingInv[i]*decayProbability[i]*supermarketContribution[i] for i in ages) + outdating*supermarketContribution[end] - holdingCosts*sum(purchasingInv[i] for i in ages) - purchasing[1]*harvestCosts[1])
	end

	@objective(model, Max, obj)
	#constraints
	if fixPurchasing
		@constraint(model, sum(purchasing[h]*harvestCosts[h]*yieldProbability[h-1] for h in yieldScenarios) >= purchasingCostContribution)
	end
	for c in wineClasses
		@constraint(model, sum(salesDiv[c,i] for i in minAgeBlending[c]:numAges) <= detDemand[c])
		@constraint(model, sum(salesDiv[c,i] for i in minAgeBlending[c]:numAges) * targetAges[c] <= sum(salesDiv[c,i] * i for i in minAgeBlending[c]:numAges))
		for i in minAgeBlending[c]:numAges
			@constraint(model, sum(usage[c,k,i] for k in blendingPatterns[c] if i in k) == salesDiv[c,i])
		end
		# @constraint(model, sum(sum(usage[c,k,i]*i for k in blendingPatterns[c] if i in k) for i in minAgeBlending[c]:numAges) == sum(salesDiv[c,i] * i for i in minAgeBlending[c]:numAges))
		for k in blendingPatterns[c]
			@constraint(model, sum(usage[c,k,i]*i for i in k) >= sum(usage[c,k,i] for i in k) * targetAges[c])
		end
	end
	@constraint(model, sum(salesDiv) <= sum(sales))
	if pU
		@constraint(model, purchasingInv[1] == sum(purchasing[h] * yieldProbability[h-1] for h in yieldScenarios))
	else
		@constraint(model, purchasingInv[1] == purchasing[1])
	end
	@constraint(model, outdating == blendingInv[ages[end]] - sales[ages[end]])
	for i in ages[1:end-1]
		@constraint(model, purchasingInv[i+1] == blendingInv[i] - sales[i])
	end
	if minAgeBlending[1] > 1
		for i in 1:minAgeBlending[1]-1
			@constraint(model, sales[i]==0)
		end
	end
	for c in 2:length(wineClasses)
		NoSalesAges = [i for i in ages if minAgeBlending[c] > i]
		for i in NoSalesAges
			@constraint(model, sales[i] == sum(salesDiv[w,i] for w in wineClasses if w < c && i >= minAgeBlending[w]))
		end
	end
	for i in minAgeBlending[wineClasses[end]]:numAges
		@constraint(model,sum(salesDiv[c,i] for c in wineClasses) == sales[i])
	end
	for i in ages
		@constraint(model, blendingInv[i] == purchasingInv[i]*(1-decayProbability[i]))
	end
	if scalePostBlending
		sF = Int(numAges / length(avgBlendingInv))
		agesBase = [i for i in 1:Int(numAges/sF)]
		for i in 1:agesBase[end]
			#ensure equal starting inventory across disaggregated age classes
			@constraint(model, blendingInv[(i-1)*sF + 1] >= avgBlendingInv[i])
			if i < agesBase[end]
				@constraint(model, purchasingInv[(i)*sF + 1] <= avgBlendingInv[i] - avgSelling[i] + 0.3)
			end
			# @constraint(model, sum(sales[j] for j in correspondenceSets[i]) == avgSelling[i])
			@constraint(model, sum(sales[j] for j in correspondenceSets[i]) >= avgSelling[i] - 0.1)
			@constraint(model, sum(sales[j] for j in correspondenceSets[i]) <= avgSelling[i] + 0.1)
		end

		@constraint(model, outdating <= 0)

		ageDeviation = ((avgSelling)' * agesBase - sum(avgServiceLevels[c]*targetAgesBase[c]*detDemand[c] for c in wineClasses))/sum(avgServiceLevels[c]*targetAgesBase[c]*detDemand[c] for c in wineClasses)

		@constraint(model, ((sales)' * ages - sum(sum(salesDiv[c,i] for i in minAgeBlending[c]:numAges)*targetAges[c] for c in wineClasses)) >= ageDeviation * sum(sum(salesDiv[c,i] for i in minAgeBlending[c]:numAges)*targetAges[c] for c in wineClasses))
		#@constraint(model, (blendingInv)' * ages >= minAvgAgeScaled * sum(blendingInv))
		for i in minAgeBlending[1]:numAges
			if any(i in targetBlendRange[c] for c in wineClasses)
				@constraint(model, sales[i] == constant + sum(factors[k] * (i^k) for k in 1:length(factors)))
				# @constraint(model, sales[i] <= scaledSelling[i]+allowedFluctuation)
				# @constraint(model, sales[i] >= scaledSelling[i]-allowedFluctuation)
			end
			for c in wineClasses
				if i >= minAgeBlending[c]
					@constraint(model, salesDiv[c,i] <= maxSales[c][i])
				end
			end
		end
		for p in yieldScenarios
			@constraint(model, purchasing[p] <= avgPricePurchasing[p-1])# + allowedFluctuation)
			@constraint(model, purchasing[p] >= avgPricePurchasing[p-1])# - allowedFluctuation)
		end
		for c in wineClasses
			if c != wineClasses[end]
				@constraint(model,sum(salesDiv[c,i] for i in minAgeBlending[c]:numAges) >= avgServiceLevels[c]*detDemand[c])
			end
			for k in blendingPatterns[c]
				#if length(k) > 2
				#	@constraint(model, sum(usage[c,k,i] for i in k) == 0)
				if length(k) == 2
					@constraint(model, usage[c,k,k[1]] == usage[c,k,k[2]])
				end
				if length(k) == 3
					@constraint(model, usage[c,k,k[1]] == usage[c,k,k[2]])
					@constraint(model, usage[c,k,k[2]] == usage[c,k,k[3]])
				end
			end
		end
	end

	optimize!(model)

	if termination_status(model) == MOI.OPTIMAL
		optimal_sales = value.(sales)
		println("sales: ", optimal_sales)
		optimal_purchasing = value.(purchasing)
		println("purchasing: ", string((optimal_purchasing)' * [yieldProbability[i-1] for i in yieldScenarios]))
		optimal_outdating = value.(outdating)
		println("outdating: ", optimal_outdating)
		optimal_salesDiv = value.(salesDiv)
		optimal_salesDiv = Dict(c=>[optimal_salesDiv[c,i] for i in minAgeBlending[c]:numAges] for c in wineClasses)
		println("salesDiv: ", optimal_salesDiv)
		optServiceLevels = [sum(optimal_salesDiv[c])/detDemand[c] for c in wineClasses]
		println("service level: ", optServiceLevels)
		optimal_purchasingInv = value.(purchasingInv)
		println("purchasingInv: ", optimal_purchasingInv)
		optimal_blendingInv = value.(blendingInv)
		optimal_objective = objective_value(model)
		println(optimal_objective)

		if scalePostBlending
			return Dict("Sales"=>optimal_sales, "Purchasing"=>optimal_purchasing, "Outdating"=>optimal_outdating, "salesDiv"=>optimal_salesDiv, "purchasingInv"=>optimal_purchasingInv, "blendingInv"=>optimal_blendingInv, "serviceLevel"=>optServiceLevels)
		else
			return Dict("Sales"=>optimal_sales, "Purchasing"=>optimal_purchasing, "salesDiv"=>optimal_salesDiv, "blendingInv"=>optimal_blendingInv, "ub" => optimal_objective, "purchasingInv"=>optimal_purchasingInv)
		end

	elseif termination_status(model) == MOI.TIME_LIMIT && has_values(model)
		suboptimal_solution = value.(sales)
		suboptimal_objective = objective_value(model)
		error("Suboptimal solution.")
	else
		error("The model was not solved correctly.")
	end
end
#----------------------------------#
optimizeAverage()
