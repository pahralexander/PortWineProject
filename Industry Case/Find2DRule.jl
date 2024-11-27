using JSON
using Distributions
using Dates
using Base.Threads
using Random
#using PyPlot
using StatsBase
using DataFrames
#using GLM
using Statistics
using XLSX
#-------------------------------------------#

#-------------------------------------------#
function find2DRule()
	for nP in [2]
		for sF in [25/6.0]
			for pU in [true]
				pSrange = pU ? [20] : ["na"]
				for sprRev in [3]
					for oDP in [0.25]
						for pS in pSrange
                            pUString = pU ? "PU" : "NoPU"
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
                            oResults = Dict()
                            mdpInstance = simulationInstance * string("/Results/",nP,"nP_",sF,"sF_",pUString,"pU_",sprRev,"sprRev_",oDP,"oDP_",pS,"pS_MDP.json")
                            open(mdpInstance, "r") do f
                                oResults = JSON.parse(f)
                            end

                            println(mdpInstance)

                            stateSpaceSize = length(keys(oResults["purchasingStateSpace"]))
                            stateSpace = Dict(s => oResults["purchasingStateSpace"][string(s)] for s in 1:stateSpaceSize)

                            optimalPurchasing = [oResults["optimalPurchasing"][string(s)] for s in 1:stateSpaceSize]
                            stateInv = [sum(stateSpace[s][2:end]) for s in 1:stateSpaceSize]
                            stateInv3 = [sum(stateSpace[s][2:3]) for s in 1:stateSpaceSize]
                            statePrice = [stateSpace[s][1] for s in 1:stateSpaceSize]
                            orderUpTo = [sum(stateSpace[s][2:end]) + optimalPurchasing[s] for s in 1:stateSpaceSize]
                            orderUpTo3 = [sum(stateSpace[s][2:3]) + optimalPurchasing[s] for s in 1:stateSpaceSize]


                            sResults = Dict()
                            simulationInstance = simulationInstance * string("/Analysis/Simulation_Aggregated/",nP,"nP_",sF,"sF_",pUString,"pU_",sprRev,"sprRev_",oDP,"oDP_",pS,"pS_Simulation.json")
                            open(simulationInstance, "r") do f
                                sResults = JSON.parse(f)
                            end
                            stateVisits = sResults["visitedPurchasing"]

                            nPriceLevels = pU ? 9 : 1
                            priceLevels = [i for i in 0:nPriceLevels-1]
                            statesWOPrice = stateSpaceSize/nPriceLevels

                            p2DRules = Dict(p=>[15,12] for p in priceLevels)
                            for p in priceLevels
                                best2D = [15,12]
                                quadFit = Inf
                                for uLevel in [i for i in 10:38]
                                    for lLevel in [i for i in 4:min(28,uLevel)]
                                        newQuadFit = 0
                                        for s in Int(1+statesWOPrice*p):Int(statesWOPrice*(p+1))
                                            newQuadFit += (max(0,uLevel-stateInv[s],lLevel-stateInv3[s])-optimalPurchasing[s])^2 * stateVisits[s]
                                        end
                                        if newQuadFit < quadFit
                                            best2D = [uLevel, lLevel]
                                            quadFit = newQuadFit
                                        end
                                    end
                                end
                                p2DRules[p] = best2D
                                println("$p done, best 2D rule: ", string(best2D))
                            end
                            sResults["2DRule"] = p2DRules
                            open(simulationInstance, "w") do f
                                JSON.print(f,sResults,4)
                            end
                        end
                    end
                end
            end
        end
    end
end
#-------------------------------------------#
find2DRule()
