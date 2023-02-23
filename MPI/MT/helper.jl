
#helper: sample Pi from Dirichlet distributions
function mysamplePi(deltaArray,BigPi)
      temp = hcat(deltaArray...) #convert array of vectors to matrix
      nLoci_array = [sum(mean(abs.(temp .- i'), dims=2) .== 0.0) + 1 for i in keys(BigPi)]
      tempPi = rand(Dirichlet(nLoci_array))
      for (i,v) in zip(keys(BigPi), tempPi)
            BigPi[i] = v
      end
      return BigPi
end


#helper: sample variances
function sample_variance(ycorr_array, nobs, df, scale, invweights, constraint)
      invweights = (invweights === false) ? false : Diagonal(invweights)
      ntraits = length(ycorr_array)
      SSE   = zeros(ntraits,ntraits) 
      for traiti in 1:ntraits
            ycorri = ycorr_array[traiti]
            for traitj in traiti:ntraits
                  ycorrj = ycorr_array[traitj]
                  SSE[traiti,traitj] = (invweights === false) ? dot(ycorri,ycorrj) : ycorri'*invweights*ycorrj
                  if constraint == true #diagonal elements only
                  break
                  end
                  SSE[traitj,traiti] = SSE[traiti,traitj]
            end
      end
      if constraint
            R = similar(SSE)
            for traiti in 1:ntraits
                  R[traiti, traiti] = (SSE[traiti, traiti] + df * scale[traiti]) / rand(Chisq(nobs + df))
            end
      else
            R  = rand(InverseWishart(df + nobs, convert(Array,Symmetric(scale + SSE))))
      end
      return R
end

