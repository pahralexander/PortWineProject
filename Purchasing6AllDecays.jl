using JSON
using Distributions
using Dates
using Base.Threads
using Random
using PyPlot
using StatsBase
HARVEST_YIELD_MEAN=4.8; OVERALL_DECAY=0.2
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
    open(string("C:/Users/ga84cib/Documents/Julia/",HARVEST_YIELD_MEAN,"hy_",OVERALL_DECAY,"od_AllDecaysAllActions.json"), "r") do f
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

	#MODEL
	#define actions
	actionSpace = [i for i in 0:nInventoryLevels-1]

	#define all possible decay transitions based on all possible states
	start = now()

	numAfterSalesStates = nInventoryLevels^(numAges-1); afterSalesStates = Dict{Int32,Array{Float64,1}}();  expectedStateReward = Dict{Int32,Float64}()
	sizehint!(expectedStateReward,stateSpaceSize); sizehint!(afterSalesStates,stateSpaceSize)
	println("Decay transition started ", now()-start)
	for s in 1:stateSpaceSize
		tempTrans = getDecays(1, s-1, stateSpaceSize, combinedDecay, nInventoryLevels, numAges)
		nTrans = size(tempTrans,1)
		expectedStateReward[s] = sum(((stateSpace[Int(s-tempTrans[i,1]+1)])' *supermarketContribution + stateContribution[Int(tempTrans[i,1])]) * tempTrans[i,2]  for i in 1:nTrans)
		afterSalesStates[s] = zeros(Float64,numAfterSalesStates)
		@inbounds @simd for i in 1:nTrans
			afterSalesStates[s][proximusState[Int(tempTrans[i,1])]] += tempTrans[i,2]
		end
		if !isapprox(sum(tempTrans[:,2]),1,atol=0.01) error("probs don't add up to 1, but ", sum(tempTrans[:,2]), " at a length of ", nTrans, " at s = ", s); end
		if s % 1000 == 0
			println(s, " transition states done in: ",now()-start)
		end
	end
	# decayTransitions = Dict();
    # open("C:/Users/ga84cib/Documents/Julia/Tmp.json", "r") do f
	# 	decayTransitions = JSON.parse(f)
	# end
	# afterSalesStates = decayTransitions["afterSalesStates"]
	# expectedStateReward = decayTransitions["expectedStateReward"]
	final = now()-start; tempTrans = nothing; reversedStateSpace = nothing; nTrans=nothing
	println("Time to create transition space: ", final)
	println(varinfo())

	# tempResults = Dict("afterSalesStates"=>afterSalesStates, "expectedStateReward"=>expectedStateReward)
	# open("C:/Users/ga84cib/Documents/Julia/Tmp.json","w") do f
	# 	JSON.print(f,tempResults,4)
	# end
	# tempResults = nothing
	# VALUE ITERATION ALGORITHM
	optimalStateActions = valueIteration(afterSalesStates,expectedStateReward,harvestCosts,actionSpace,stateSpaceSize,yieldProbability,nInventoryLevels)

	resultsDict = Dict("stateSpace"=>stateSpace,"actionSpace"=>actionSpace,"expectedStateReward"=>expectedStateReward,"optimalStateActions"=>optimalStateActions)
	open("C:/Users/ga84cib/Documents/Julia/ResultsPurchasing6AllDecays_",HARVEST_YIELD_MEAN,"hy_",OVERALL_DECAY,"od.json","w") do f
		JSON.print(f,resultsDict,4)
	end
end

#----------------------------------------#
# value iteration algorithm
function valueIteration(afterSalesStates,expectedStateReward,harvestCosts,actionSpace,stateSpaceSize,yieldProbability,nInventoryLevels)
	priceIdentifier = Int(stateSpaceSize/nInventoryLevels)
	yieldSections = [i*priceIdentifier+1:i*priceIdentifier+priceIdentifier for i in 0:nInventoryLevels-1]
	preActionState = [(s-1)%priceIdentifier+1 for s in 1:stateSpaceSize]
	statePrice = [harvestCosts[Int(floor(s/priceIdentifier))+1] for s in 0:stateSpaceSize-1]
	start=now()
	convergence = Inf
	epsilon = 0.001
	Vs = zeros(Float64,2,stateSpaceSize)
	policy = zeros(Int32,stateSpaceSize)
	iteration = 0
	while convergence > epsilon
		iteration += 1
		println("Iteration: ", iteration)
		previousValues = [view(Vs,(iteration-1)%2+1,yieldSections[y]) * yieldProbability[y] for y in 1:nInventoryLevels]
		expectedPreviousValues = Dict(); sizehint!(expectedPreviousValues, stateSpaceSize)
		for i in 1:stateSpaceSize
			expectedPreviousValues[i] = sum((afterSalesStates[i])' * previousValues[y] for y in 1:nInventoryLevels)
		end
		@inbounds for s in Int32(1):Int32(stateSpaceSize)
			Vs[iteration%2+1,s] = -Inf
			for a in actionSpace
				afterActionState = a*priceIdentifier + preActionState[s]
				# transValues = view(Vs,iteration%2+1,tStates[transType[s][a]])
				#Value=
				ActionValue = -(statePrice[s]*a) + expectedStateReward[afterActionState] + expectedPreviousValues[afterActionState]
				if ActionValue > Vs[iteration%2+1,s]
					Vs[iteration%2+1,s] = ActionValue
					policy[s] = a
				end
				if Vs[iteration%2+1,s] == -Inf
					error("Value not updated: ", ActionValue, " at actionContribution: ", i[1], " for action: ", a, " ", ActionValue>Vs[iteration%2+1,s])
				end
			end
			if s%20000 == 0
				println(s, " done in ", now() - start, "Value: ", Vs[iteration%2+1,s])
			end
		end
		M = maximum(view(Vs,iteration%2+1,:) - view(Vs,(iteration-1)%2+1,:))
		m = minimum(view(Vs,iteration%2+1,:) - view(Vs,(iteration-1)%2+1,:))
		println(M)
		println(m)
		convergence = M-m
		println("convergence: ",convergence)
		println("Time: ",now()-start)
	end
	final = now() - start
	println("Solved in ", final)
	println("Completed in ", iteration, " iterations")

	return Dict(i => (policy[i],Vs[iteration%2+1,i]) for i in 1:stateSpaceSize)

end
#----------------------------------------#
function getDecays(seed,state,posMultiplier,combinedDecay,nInventoryLevels,numAges,transState=1,transProbability=1.0)
	elementDivisor = copy(posMultiplier)/nInventoryLevels
	currentElement = copy(seed)
	if currentElement >= numAges+1
		return [transState transProbability]
	end
	transitionSpace = Array{Float64}(undef,0,2)
	numInStock = Int8(floor(state/elementDivisor))
	if numInStock > 0
		for dec in 0:numInStock
			newTransState = transState + copy(elementDivisor)*(numInStock-dec)
			newTransProb = transProbability * combinedDecay[numInStock][currentElement,dec+1]
			newTransitions = getDecays(seed+1,state%elementDivisor,elementDivisor,combinedDecay,nInventoryLevels,numAges,newTransState,newTransProb)
			transitionSpace = vcat(transitionSpace,newTransitions)
		end
	else
		newTransitions = getDecays(seed+1,state%elementDivisor,elementDivisor,combinedDecay,nInventoryLevels,numAges,transState,transProbability)
		transitionSpace = vcat(transitionSpace,newTransitions)
	end
	return transitionSpace
end
#------------------------------------------------#
main()
