function TrackingPlot(TckResult,prnList)
    for i=1:length(prnList)
        prnNumber = prnList(i);
        figure(prnNumber);

        subplot(4,3,1);
        plot(TckResult(prnNumber).P_i(:),TckResult(prnNumber).P_q(:),'b*');
        title('Discrete-Time Scatter Plot');
        xlabel('I prompt');
        ylabel('Q prompt');
        grid on;

        subplot(4,3,2:3);
        plot(TckResult(prnNumber).P_i(:),'r-');
        hold on;
        plot(TckResult(prnNumber).P_q(:),'b-');
        title('Ip and Qp of Tracking Result');
        legend('I_P','Q_P');
        grid on;

        subplot(4,3,4);
        plot(TckResult(prnNumber).PLLdiscri(:),'b-');
        title('PLL Discriminator');
        xlabel('Time(ms)');
        ylabel('Amplitude');
        grid on;

        subplot(4,3,5:6)
        plot(sqrt(TckResult(prnNumber).E_i(:).^2+TckResult(prnNumber).E_q(:).^2),'b*-');hold on
        plot(sqrt(TckResult(prnNumber).L_i(:).^2+TckResult(prnNumber).L_q(:).^2),'g*-');hold on
        plot(sqrt(TckResult(prnNumber).P_i(:).^2+TckResult(prnNumber).P_q(:).^2),'r*-');hold on
        legend('Early','Late','Prompt');
        title('Correlation');
        xlabel('Time(ms)');
        grid on;

        subplot(4,3,7);
        plot(TckResult(prnNumber).carrFreq,'b');
        title('Carrier frequency');
        xlabel('Time(ms)');
        ylabel('Hz');
        grid on;

        subplot(4,3,8);
        plot(TckResult(prnNumber).DLLdiscri,'r-');
        title('DLL Discriminator');
        xlabel('Time(ms)');
        ylabel('Amplitude');
        grid on;

        subplot(4,3,9);
        plot(TckResult(prnNumber).codedelay,'r*-');
        title('Code Delay Sample');
        xlabel('Time(ms)');
        ylabel('Sample');
        grid on;

        RawNavigationData = TckResult(prnNumber).P_i(:);
        NaviData(find(RawNavigationData>=0)) = 1;
        NaviData(find(RawNavigationData<0)) = -1; 
        datalength = length(TckResult(prnNumber).P_i);
        subplot(4,3,10:12);
        stairs(NaviData);
        axis([1 datalength -1.2 1.2])
        title('Navigation Data');
        xlabel('Time(ms)');
        grid on;    
    end