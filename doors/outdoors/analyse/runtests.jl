n = 100
xy = DimVector([SV(0, i) for i in 1:n], Ti(range(2, 50, n)))
intervention = 10
xy[Ti = intervention..Inf] .+= Ref(SV(4,4))

glue_intervention!(xy, intervention)
@test all(<(1.001) ∘ norm, diff(xy))
