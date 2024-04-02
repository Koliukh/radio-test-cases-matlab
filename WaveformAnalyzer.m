classdef WaveformAnalyzer < handle
    %% Описание класса
    %
    % 1. Класс читает данные (во временной области) на выходе OFDM модулятора сигнала, а также информацию о параметрах формирователя
    %
    % 2. Строит метрики: спектральная плотность мощности в частотной области, графическое представление созвездия на комплексной плоскости,
    % среднеквадратичное значение модуля вектора ошибки (EVM)
    %
    % Входные данные:
    %
    % waveformSource - массив содержащий отчеты baseband сигнала во временной области на выходе OFDM модулятора
    %
    % waveformInfo - структура с параметрами OFDM модулятора и пейлоуда:
    %       Nfft               - кол-во спектрально-временных отчетов дискретного преобразования Фурье
    %       SampleRate         - частота семплирования [Гц]
    %       CyclicPrefixLengths/SymbolLengths - длины циклического преффикса и OFDM символов [кол-во временных отчетов]
    %       SymbolsCount       - кол-во символов на слот радиокадра
    %       subCarriersCount   - кол-во поднесущих
    %       payloadSymbols     - информационные символы
    %       payloadSymbolsIdxs - индексы ресурсных элементов отведенные для передачи payloadSymbols
    %
    % Поля класса:
    %
    %       rmsEvm            - среднеквадратичное значение модуля вектора ошибки
    %       waveformMeanPower - среднеквадратичное значение мощности сигнала
    %       channelBandwidth  - ширина полосы канала
    %       noiseMeanPower    - среднеквадратичное значение мощности шума
    %       modulationType    - тип модуляционной схемы
    %       waveformDuration  - длина анализируемого сигнала
    %

    properties
        Nfft
        SampleRate
        CyclicPrefixLengths
        SymbolLengths
        SymbolsCount
        Windowing
        subCarriersCount
        SymbolPhases
        SymbolsPerSlot
        symbolsCount
        payloadSymbols
        payloadSymbolsIdxs
        rxWaveform
        rx_sym
    end
    
    properties
        rmsEvm
        waveformMeanPower
        channelBandwidth
        noiseMeanPower
        modulationType
        waveformDuration
        dopplershift
    end

    methods
        function this = WaveformAnalyzer()
            % Конструктор класса. Чтение waveform-ы во временной области и структуры с информацией
            % необходимой для дальнейшей обработки данных и заполнения полей класса
            load('waveformInfo.mat')
            load('waveformSource.mat')
            
            this.Nfft = info.Nfft;
            this.SampleRate = info.SampleRate;
            this.CyclicPrefixLengths = info.CyclicPrefixLengths;
            this.SymbolLengths = info.SymbolLengths;
            this.symbolsCount = info.symbolsCount;
            this.Windowing = info.Windowing ;
            this.SymbolPhases = info.SymbolPhases ;
            this.SymbolsPerSlot = info.SymbolsPerSlot ;
            this.symbolsCount = info.symbolsCount ;
            this.subCarriersCount = info.subCarriersCount;
            this.payloadSymbols = info.payloadSymbols;
            this.payloadSymbolsIdxs = info.payloadSymbolsIdxs;
            this.rxWaveform = rxWaveform;
        end

        function calcWaveformParameters(this)
            tx = this.payloadSymbols;
            N = this.Nfft - this.subCarriersCount;
            ofdmDemod = comm.OFDMDemodulator('FFTLength',this.Nfft,...
                'NumGuardBandCarriers',[N/2;N/2],...
                'CyclicPrefixLength',this.CyclicPrefixLengths,...
                'NumSymbols',this.symbolsCount,'NumReceiveAntennas',1);
            % disp(ofdmDemod)
            % showResourceMapping(ofdmDemod)
            [dataOut] = step(ofdmDemod,this.rxWaveform).';
            rx_ = dataOut(:);
            rx = rx_(this.payloadSymbolsIdxs);
            tx = this.payloadSymbols;
            mem=500; 
            % тут видимо нужно как то правильно канал учесть, потому что мой способ очевидно не правильный
            rx_eq = this.lin_eq(rx, tx, mem);
            % figure; plot(xcorr(real(tx),real(rx_eq)))
            % figure; hold on; plot(rx,'.') ;plot(tx,'.');plot(tx-rx,'.');
            this.rx_sym = rx_eq;
            this.waveformMeanPower = rms(this.rxWaveform);
            this.channelBandwidth
            this.noiseMeanPower = rms(tx - rx_eq);
            this.modulationType  = 'QAM64';
            this.waveformDuration  = length(this.rxWaveform);
        end

        function calcdopplerSHift

        end

        function plotPowerSpectrumDensity(this)
            [Pxx,f] = pwelch(this.rxWaveform,[],[],[],this.SampleRate,'centered');
            figure; hold on; grid on;
            plot(f/1e6,10*log10(Pxx));
            xlabel('Frequency (MHz)')
            ylabel('Power (dB)')
            legend('pwelch')
            title(['Fs=', mat2str(this.SampleRate/1e6),' MHz'])
        end

        function plotPayloadConstellation(this)
            figure; hold on; grid on; title('constellation')
            plot(this.payloadSymbols,'r*');
            xlabel('Re'); ylabel('Im');
        end

        function calcEvmPerformance(this)
           refSym = this.payloadSymbols;
            rxSym = this.rx_sym; 
            const  = unique(refSym);
            len = length(refSym);
            e_re = real(refSym-rxSym);
            e_im = imag(refSym-rxSym);
            sig_pwr = sum(real(refSym).^2 +imag(refSym).^2)/len;
            aver_const_pwr = mean(abs(const).^2); 
            peak_const_pwr = max(abs(const).^2);

            EVMref = 100*sqrt(sum(e_re.^2 + e_im.^2)/len/sig_pwr);
            EVMaver = 100*sqrt(sum(e_re.^2 + e_im.^2)/len/aver_const_pwr);
            EVMpeak = 100*sqrt(sum(e_re.^2 + e_im.^2)/len/peak_const_pwr);
            fprintf('EVMref = %1.3f  EVMaver = %1.3f  EVMpeak = %1.3f \n',EVMref,EVMaver,EVMpeak);
            this.rmsEvm = [EVMref,EVMaver,EVMpeak];
            %NMSE
            %nmse = 10*log10(sum(abs(tx-rx_eq).^2)/length(tx))
        end
        function eq_out = lin_eq(this,rx, tx, mem)
            rx = rx(:);
            tx = tx(:);
            V_= [];
            for shift = -mem:mem
                V_ = [V_,circshift(rx,shift)];
            end
            V = [real(V_),imag(V_),abs(V_)];
            coeff = V\tx;
%             figure; stem(coeff)
            eq_out = V*coeff;
        end
    end
end
