function SR_SYSTEM = IEEELetters2017_System_B()


%% Create soundfield reproduction system structure
SR_SYSTEM = Current_Systems.IEEELetters2017_System_A;
SR_SYSTEM.Main_Setup(3:end)  = [];

end
