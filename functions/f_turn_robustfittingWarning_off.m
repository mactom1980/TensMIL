function [] = f_turn_robustfittingWarning_off(x)

if x
    warning('off','stats:statrobustfit:IterationLimit')
else
    warning('on','stats:statrobustfit:IterationLimit')
end
end