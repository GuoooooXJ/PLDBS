function [stLabjack,ljudObj,ljhandle,channel]=Labjack_Ini()
    % Make the UD .NET assembly visible in MATLAB.
    ljasm = NET.addAssembly('LJUDDotNet');
    ljudObj = LabJack.LabJackUD.LJUD;
    
    % Read and display the UD version.
    disp(['UD Driver Version = ' num2str(ljudObj.GetDriverVersion())])
    
    try 
        % Open the first found LabJack U3.
        [ljerror, ljhandle] = ljudObj.OpenLabJackS('LJ_dtU3', 'LJ_ctUSB', '0', true, 0);    
    catch e
        if(isa(e, 'NET.NetException'))
            eNet = e.ExceptionObject;
            if(isa(eNet, 'LabJack.LabJackUD.LabJackUDException'))
                disp(['UD Error: ' char(eNet.ToString())])
            else
            disp(['.NET Error: ' char(eNet.ToString())])
            end
        end
        disp(getReport(e))
        stLabjack=-1;
        return;
    end
    
    % Set DAC0 to vx volts.    
    % Set the state of FIO7.
    channel = 7; 
    state = 0; 
    ljudObj.eDO(ljhandle, channel, state);
    
    stLabjack=1;
    
