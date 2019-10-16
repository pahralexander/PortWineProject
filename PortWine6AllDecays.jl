# imports
using JSON
using Distributions
using Dates
using Base.Threads
using Random
using PyPlot
using StatsBase
@show Threads.nthreads()
#-------------------------------------------#
function main()
	start = now()
	#initialize Age classes
	ages = [1,2,3,4,5,6]
	numAges = ages[end]

	#initialize Wine classes
	wineClasses = ["young", "old"]
	wineClassIndex = Dict("young"=>1,"old"=>2)
	targetAges = Dict("young"=>3, "old"=>5)
	detDemand = Dict("young"=>2, "old"=>2)
	totalDemand = sum(values(detDemand))
	#initialize Inventory levels
	nInventoryLevels = 7

	# define maximum sales for each age class
	maxSales = Dict(i=>0 for i in ages)
	for i in ages
		if i >= targetAges["old"]
			maxSales[i] = min(totalDemand,nInventoryLevels-1)
		else
			for j in (detDemand["old"]-1):-1:1
				if j*i+(detDemand["old"]-j)*ages[end] >= targetAges["old"]*detDemand["old"]
					maxSales[i] += j
					break
				end
			end
			if i >= targetAges["young"]
				maxSales[i] += min(detDemand["young"],nInventoryLevels-1)
			else
				for j in (detDemand["young"]-1):-1:1
					if j*i+(detDemand["young"]-j)*ages[end] >= targetAges["young"]*detDemand["young"]
						maxSales[i] += j
						break
					end
				end
			end
		end
	end
	# println("maxSales ",maxSales)

	#define minimum inventory levels in each age class if only degradation of the whole class is considered
	minYield = 3
	#define minimum sum of all ages for different demand scenarios
	minimumAgeSum=Dict(i=>0 for i in 0:totalDemand)
	for i in keys(minimumAgeSum)
		remYoung=detDemand["young"]
		remOld=detDemand["old"]
		c=i
		while c>0
			c-=1
			if remYoung>0
				minimumAgeSum[i]+=targetAges["young"]
				remYoung-=1
			else
				minimumAgeSum[i]+=targetAges["old"]
			end
		end
	end
	#define possible age patterns for different sales quantities
	possiblePatterns = Dict()
	for i in 0:totalDemand
		possiblePatterns[i]=[]
		for k in 0:(min(detDemand["young"],i))
			if (i-k) <= detDemand["old"]
				append!(possiblePatterns[i],[(k,i-k)])
			end
		end
	end
	println("possiblePatterns ",possiblePatterns)

	#DATA(key P&L figures)
	#price of wine in respective age classes
	sellingPrice = Dict("young"=>500, "old"=> 5000)
	supermarketPrice = [200+(i-1)*20 for i in 1:ages[end]]
	#bottling costs of each respective age class
	bottlingCosts = Dict("young"=>90, "old"=>1200)
	supermarketContribution = supermarketPrice .- bottlingCosts["young"]
	#wine costs (=production costs) of each respective age class
	# wineCosts = Dict("young"=> 200, "old"=>273)
	#holding costs
	holdingCosts = 80
	overallDecay = 0.2

	#define harvest yield probabilities (truncating a discrete normal distribution between 3 and 6)
	muH = 4.8
	sigmaH = 1.2
	truncDist = Truncated(Normal(muH,sigmaH),minYield-0.5,nInventoryLevels-0.5)
	yieldProbability = Dict(i=>cdf(truncDist,i+0.5)-cdf(truncDist,i-0.5) for i in minYield:nInventoryLevels-1)
	# println(yieldProbability)
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
	#define states
	start = now()
	allStates = defineStateSpace(0,numAges,nInventoryLevels)
	stateSpaceSize = length(allStates[:,1])
	stateSpace = Dict{Int32,Array{Int8,1}}(); sizehint!(stateSpace,stateSpaceSize)
	stateSpace = Dict(s=>allStates[s,:] for s in 1:stateSpaceSize);
	reversedStateSpace = Dict(v=>k for (k,v) in stateSpace)
	println("Number of states: ", stateSpaceSize)
	final=now()-start
	println("Time to create state space: ", final)

	#define all possible transitions based on all possible states
	start = now()
	transInputs = collect(Set(values(Dict(i=>allStates[i,1:end-1] for i in 1:stateSpaceSize)))); allStates = nothing;
	tProbs = Dict{Array{Int8,1},Array{Float64,1}}(); tSupermarket = Dict{Array{Int8,1},Float64}()
	sizehint!(tProbs,nInventoryLevels^(numAges-1)); sizehint!(tSupermarket,nInventoryLevels^(numAges-1))
	println("Transition generation started after ", now()-start)
	for t in 1:length(transInputs)
		input = copy(transInputs[t])
		tempTrans = getTransitions(Int8(1),input,numAges,combinedDecay,yieldProbability)
		nTrans = size(tempTrans,1)
		tStates = [reversedStateSpace[tempTrans[i,1:numAges]] for i in 1:nTrans]
		# tProbs[transInputs[t]] = tempTrans[:,numAges+1]
		tDecay = view(tempTrans,:,numAges+2:numAges*2+1) * supermarketContribution
		tSupermarket[input] = ((view(tempTrans,:,numAges+1))' * tDecay)[1]
		tProbs[input] = zeros(Float64,stateSpaceSize)
		@inbounds @simd for s in 1:nTrans
			tProbs[input][tStates[s]] += tempTrans[s,numAges+1]
		end
		if !isapprox(sum(tProbs[transInputs[t]]),1,atol=0.001) error("probs don't add up to 1, but ", sum(tProbs[transInputs[t]]), " at a length of ", length(tProbs[transInputs[t]]), " at t = ", t); end
		if t % 1000 == 0
			println(t, " transition states done in: ",now()-start)
		end
	end
	final = now()-start; tempTrans = nothing; transInputs = nothing; reversedStateSpace = nothing; tStates = nothing; nTrans=nothing
	println("Time to create transition space: ", final)
	println(varinfo())

	#defineActions
	start=now()
	allActions = defineActionSpace(1,maxSales,numAges,totalDemand,minimumAgeSum,possiblePatterns,targetAges,detDemand)
	actionSpaceSize = length(allActions[:,1])
	actionSpace = Dict(a=>allActions[a,1:numAges] for a in 1:actionSpaceSize)
	demandFulfillment = Dict(a=>allActions[a,numAges+1:end] for a in 1:actionSpaceSize)
	fulfillYoungOld = Dict(a => Dict("young"=>demandFulfillment[a][1], "old"=>demandFulfillment[a][2]) for a in 1:actionSpaceSize)
	println("Number of actions: ", actionSpaceSize)
	actionContribution = Dict(a=>sum(fulfillYoungOld[a][i]*(sellingPrice[i]-bottlingCosts[i]) for i in wineClasses) for a in 1:actionSpaceSize); allActions=nothing
	#define applicable actions for each state
	stateActions  = Dict{Int32,Array{Int16,1}}()
	sizehint!(stateActions,stateSpaceSize)
	for s in Int32(1):Int32(stateSpaceSize)
		stateActions[s] = Array{Int16,1}()
		stateHighestDemand = 0
		for a in Int16(1):Int16(actionSpaceSize)
			totalAction = sum(demandFulfillment[a])
			if isCompatible(stateSpace[s],actionSpace[a]) && totalAction >= stateHighestDemand
				#actions may only be executed if demand is satisfied to the greatest extent possible
				# if totalAction > stateHighestDemand
				# 	stateActions[s] = Array{Int16,1}()
				# 	stateHighestDemand = totalAction
				# end
				append!(stateActions[s],a)
			end
		end
	end
	countActions = 0
	for s in 1:stateSpaceSize
		countActions += length(stateActions[s])
	end
	println("Avg number of actions per state: ", countActions/stateSpaceSize)
	final = now() - start; demandFulfillment=nothing
	println("Time to create state specific actions space: ", final)

	#define immediate reward per state-action pair
	start = now()
	stateActionInfo = getRewardAndTransType(stateSpaceSize,actionContribution,stateActions,tSupermarket,holdingCosts,supermarketContribution[end],stateSpace,actionSpace)
	final = now() - start
	println("stateActionInfo: ", sizeof(stateActionInfo), " transProbs: ", sizeof(tProbs)); tSupermarket = nothing;
	println("Time to create state/action specific transitions: ",final)

	# VALUE ITERATION ALGORITHM
	optimalStateActions,numIterations = valueIteration(stateActionInfo,tProbs,stateSpaceSize)
	optimalPolicyRewards = Dict(s=>stateActionInfo[s][optimalStateActions[s][1]][1] for s in 1:stateSpaceSize)

	resultsDict = Dict("stateSpace"=>stateSpace,"actionSpace"=>actionSpace,"fulfillYoungOld"=>fulfillYoungOld,"stateActions"=>stateActions,"optimalStateActions"=>optimalStateActions, "optimalPolicyRewards"=>optimalPolicyRewards, "actionContribution"=>actionContribution, "numIterations"=>numIterations)
	open(string("C:/Users/ga84cib/Documents/Julia/",muH,"hy_",overallDecay,"od_AllDecaysAllActions.json"),"w") do f
		JSON.print(f,resultsDict,4)
	end
end

#----------------------------------------#
# value iteration algorithm
function valueIteration(stateActionInfo,transProbs,stateSpaceSize)
	start=now()
	convergence = Inf
	epsilon = 0.001
	Vs = zeros(Float64,2,stateSpaceSize)
	policy = zeros(Int32,stateSpaceSize)
	iteration = 0
	while convergence > epsilon
		iteration += 1
		println("Iteration: ", iteration)
		previousValues = view(Vs,(iteration-1)%2+1,:)
		expectedPreviousValues = Dict(); sizehint!(expectedPreviousValues,length(keys(transProbs)))
		for t in keys(transProbs)
			expectedPreviousValues[t] = (transProbs[t])' * previousValues
		end
		@inbounds for s in Int32(1):Int32(stateSpaceSize)
			Vs[iteration%2+1,s] = -Inf
			for (a,i) in stateActionInfo[s]
				# transValues = view(Vs,iteration%2+1,tStates[transType[s][a]])
				#Value=
				ActionValue = i[1] + expectedPreviousValues[i[2]]
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

	return Dict(i => (policy[i],Vs[iteration%2+1,i]) for i in 1:stateSpaceSize),iteration

end
#----------------------------------------#
#get the reward and the transition type of a certain state&action pair
function getRewardAndTransType(stateSpaceSize,actionContribution,stateActions,tSupermarket,holdingCosts,outdatePrice,stateSpace,actionSpace)
	start = now()
	println("Transitions done for ")
	stateActionInfo = Dict{Int32,Dict{}}()
	sizehint!(stateActionInfo,stateSpaceSize)
	@inbounds for s in Int32(1):Int32(stateSpaceSize)
		stateActionInfo[s] = Dict{Int16,Tuple{Float64,Array{Int8,1}}}()
		lsa = length(stateActions[s])
		sizehint!(stateActionInfo[s],lsa)
		@inbounds for a in stateActions[s]
			#get deterministic outcomes from state-action-pair
			newState = stateSpace[s] - actionSpace[a]
			carryOnState = newState[1:end-1]
			#fetch transition information based on interim inventory state
			stateActionInfo[s][a] = actionContribution[a] + newState[end]*outdatePrice - sum(carryOnState)*holdingCosts + tSupermarket[carryOnState],carryOnState
		end
		if s%5000 == 0
			println(s," states in ", now()-start)
		end
	end
	return stateActionInfo
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
function defineActionSpace(seed,maxSales,numAges,totalDemand,minimumAgeSum,possiblePatterns,targetAges,detDemand,actionVector=Array{Int8,1}(zeros(numAges+length(keys(targetAges)))))
	pos = copy(seed)
	if pos >= numAges+1
		check = isReasonable(actionVector,numAges,minimumAgeSum,possiblePatterns,targetAges,detDemand)
		if check != nothing
			aV = copy(actionVector)
			for p in 1:length(check)
				aV[numAges+p] = check[p]
			end
			return transpose(aV)
		else
			return nothing
		end
	end
	actionSpace = Array{Int8}(undef,0,numAges+length(keys(targetAges)))
	currentElement = pos
	salesMax = maxSales[currentElement]
	for q in Int8(0):Int8(min(totalDemand,salesMax))
		aV = copy(actionVector)
		aV[pos] = q
		newActions = defineActionSpace(pos+1,maxSales,numAges,totalDemand-q,minimumAgeSum,possiblePatterns,targetAges,detDemand,aV)
		if newActions != nothing
			actionSpace = vcat(actionSpace,newActions)
		end
	end
	return actionSpace
end
#------------------------------------------------#
#checks whether a particular action is valid
function isReasonable(actionVector,numAges,minimumAgeSum,possiblePatterns,targetAges,detDemand)
	daV = copy(actionVector)
	salesAmount = sum(daV)
	ageSum=sum(i*daV[i] for i in 1:numAges)
	if ageSum < minimumAgeSum[salesAmount]
		return nothing
	end
	salesAges = []
	counter = 1
	#determine the ages of the wines sold in the respective action
	while counter <= numAges
		if daV[counter] > 0
			append!(salesAges,counter)
			daV[counter] -= 1
		end
		if daV[counter] <= 0
			counter += 1
		end
	end
	#check if the wine ages are compatible with a possible sales pattern
	for p in possiblePatterns[salesAmount]
		pSalesAges = copy(salesAges)
		YDem = p[1]
		FF = 0
		ageFF = 0
		counter = 0
		while YDem > 0
			if pSalesAges == [] || counter > p[1]
				return nothing
			end
			counter += 1
			for a in sort(pSalesAges)
				if a+ageFF >= targetAges["young"]*(FF+1)
					deleteat!(pSalesAges,findfirst(isequal(a),pSalesAges))
					FF += 1; ageFF += a; YDem -= 1; break
				end
			end
		end
		agesForOld = isempty(pSalesAges) ? 0 : sum(pSalesAges)
		if agesForOld >= targetAges["old"]*p[2]
			return p
		end
	end
	return nothing
end
#------------------------------------------------#
function isCompatible(state,action)
	if any(state[k] < action[k] for k in 1:length(action))
		return false
	elseif sum(state)*(1/3) < sum(action)
		return false
	else
		return true
	end
end
#------------------------------------------------#
function getTransitions(seed,newState,numAges,combinedDecay,yieldProbability,transState=vcat(zeros(numAges),[1.0],zeros(numAges)))
	currState = copy(transState)
	currentElement = copy(seed)
	if currentElement >= numAges+1
		return transpose(currState)
	end
	transitionSpace = Array{Float64}(undef,0,numAges*2+1)
	if currentElement == 1
		for y in keys(yieldProbability)
			for d in 0:y
				transState = copy(currState)
				transState[currentElement] = y-d; transState[numAges+1+currentElement]=d; transState[numAges+1] = currState[numAges+1]*yieldProbability[y]*combinedDecay[y][currentElement,d+1]
				newTransitions = getTransitions(currentElement+1,newState,numAges,combinedDecay,yieldProbability,transState)
				transitionSpace = vcat(transitionSpace,newTransitions)
			end
		end
	else
		prevInv = newState[currentElement-1]
		if prevInv > 0
			for d in 0:prevInv
				transState = copy(currState)
				transState[currentElement] = prevInv - d; transState[numAges+1+currentElement]=d; transState[numAges+1] = transState[numAges+1]*combinedDecay[prevInv][currentElement,d+1]
				newTransitions = getTransitions(currentElement+1,newState,numAges,combinedDecay,yieldProbability,transState)
				transitionSpace = vcat(transitionSpace,newTransitions)
			end
		else
			transState = copy(currState)
			newTransitions = getTransitions(currentElement+1,newState,numAges,combinedDecay,yieldProbability,transState)
			transitionSpace = vcat(transitionSpace,newTransitions)
		end
	end
	return transitionSpace
end
#------------------------------------------------#
main()
