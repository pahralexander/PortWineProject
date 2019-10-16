using JSON
using Distributions
using Dates
using Base.Threads
using Random
using PyPlot
using StatsBase
HARVEST_YIELD_MEAN = 4.5; OVERALL_DECAY = 0.1;
#-------------------------------------------#
function main()
    #initialize Age classes
	ages = [1,2,3,4,5,6]
	numAges = ages[end]

	#initialize Inventory levels
	nInventoryLevels = 7

	#initialize Wine classes
	wineClasses = ["young", "old"]
	wineClassIndex = Dict("young"=>1,"old"=>2)
	targetAges = Dict("young"=>3, "old"=>5)
	detDemand = Dict("young"=>2, "old"=>2)
	totalDemand = sum(values(detDemand))
	stateMultipliers = [nInventoryLevels^(numAges-i) for i in 1:numAges]

	start = now()

	#define minimum inventory levels in each age class if only degradation of the whole class is considered
	minYield = 0

	#DATA(key P&L figures)
	#price of wine in respective age classes
	sellingPrice = Dict("young"=>500, "old"=> 5000)
	supermarketPrice = [200+(i-1)*20 for i in 1:ages[end]]
	#bottling costs of each respective age class
	bottlingCosts = Dict("young"=>90, "old"=>1200)
	supermarketContribution = supermarketPrice .- bottlingCosts["young"]
	#wine costs (=purchasing costs) of each respective age class
	wineCostStart = 140
	costStep = 20
	harvestCosts = [wineCostStart+costStep*i for i in 0:nInventoryLevels-1]
	#holding costs
	holdingCosts = 80
	overallDecay = OVERALL_DECAY

	#load results from blending model
	results = Dict();
    open("C:/Users/ga84cib/Documents/Julia/45hy_01od_AllDecaysAllActions.json", "r") do f
		results = JSON.parse(f)
	end
	stateSpaceSize = length(values(results["stateSpace"]))
	stateSpace = Dict(parse(Int32,i)=>results["stateSpace"][i] for i in keys(results["stateSpace"]))
	salesSpace = Dict(parse(Int16,i)=>results["actionSpace"][i] for i in keys(results["actionSpace"]))
	salesSpaceSize = length(values(salesSpace))
	#encode sales vectors to state numbers (WITHOUT +1)
	stEncSalesSpace = Dict(a=>(salesSpace[a])' * stateMultipliers for a in 1:salesSpaceSize)
	fulfillYoungOld = Dict(parse(Int16,i)=>results["fulfillYoungOld"][i] for i in keys(results["fulfillYoungOld"]))
	optimalStateSales = Dict(parse(Int32,s)=>results["optimalStateActions"][s][1] for s in keys(results["optimalStateActions"]))
	stEncOptimalStateSales = [stEncSalesSpace[optimalStateSales[s]] for s in 1:stateSpaceSize]
	stateOutdating = [stateSpace[s][end] - salesSpace[optimalStateSales[s]][end] for s in 1:stateSpaceSize]
	outdatingContribution = stateOutdating .* supermarketContribution[end]
	salesContribution = [sum(fulfillYoungOld[a][i]*(sellingPrice[i]-bottlingCosts[i]) for i in wineClasses) for a in 1:salesSpaceSize]
	stateSalesContribution = [salesContribution[optimalStateSales[s]] for s in 1:stateSpaceSize]
	stateContribution = [stateSalesContribution[s] + outdatingContribution[s] for s in 1:stateSpaceSize]
	reversedStateSpace = Dict(v=>k for (k,v) in stateSpace)
	proximusState = [(stateSpace[s][1:end-1] - salesSpace[optimalStateSales[s]][1:end-1])' * stateMultipliers[2:end] + 1 for s in 1:stateSpaceSize]

	println("Number of states: ", stateSpaceSize)
	#define harvest yield probabilities for price definition (truncating a discrete normal distribution between 0 and 6)

	open("C:/Users/ga84cib/Documents/Julia/",HARVEST_YIELD_MEAN,"hy_",OVERALL_DECAY,"od_AllDecaysAllActions.json", "r") do f

	end

	muH = 3.0
	sigmaH = 1.2
	truncDist = Truncated(Normal(muH,sigmaH),minYield-0.5,nInventoryLevels-0.5)
	yieldProbability = [cdf(truncDist,i+0.5)-cdf(truncDist,i-0.5) for i in minYield:nInventoryLevels-1]
	println(yieldProbability)
	# items = [i for i in minYield:nInventoryLevels-1]
	# weights = [yieldProbability[i] for i in minYield:nInventoryLevels-1]
	# sampleYields = sample(items,Weights(weights),50000)
	# hist(sampleYields,bins=nInventoryLevels-minYield)
	# show()

	#define decay probabilities (use discrete Weibull distribution as in Nakagawa and Osaki (1975)) + Bernoulli-process resulting partial decay probabilities
	q = 0.66
	β = 0.8
	decayProbability = zeros(numAges)
	approxdecayProb = zeros(numAges)
	for k in 1:ages[end]
	    approxdecayProb[k] = q^(k^(β))-q^((k+1)^(β))
	end
	for k in 1:ages[end]
	   decayProbability[k] = round(approxdecayProb[k]/sum(values(approxdecayProb))*overallDecay,digits=3)
	end
	println(decayProbability)

	println("DECAY PROBABAILITY:\n",decayProbability, sum(values(decayProbability)))
	combinedDecay = Dict(s=>zeros(numAges,s+1) for s in 1:nInventoryLevels-1)
	for s in 1:nInventoryLevels-1
		for d in 0:s
			for a in 1:numAges
				combinedDecay[s][a,d+1] = binomial(s,d)*(decayProbability[a]^d)*((1-decayProbability[a])^(s-d))
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

	#USE BINOMIAL FOR SIMULATION RUNS

end
#-------------------------------------------#
main()
