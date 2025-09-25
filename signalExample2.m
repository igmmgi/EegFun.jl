function signalExample2
% signalExample2

% background colour and font
mygui.defaultBackground = get(0,'defaultUicontrolBackgroundColor');
mygui.defaultFixedFont  = get(0,'FixedWidthFontName');

% figure window size (open position x y height width)
mygui.size = [0 0 0.25 0.5];

% figure settings
hFigure = figure('Units', 'normalized', ...
  'outerposition', mygui.size, ...
  'Color', mygui.defaultBackground, ...
  'Visible', 'off', ...
  'Name', 'Signal Example 2', ...
  'numberTitle', 'off');

% Axes for plotting the selected plot
hPlotAxes1 = axes('Parent', hFigure, ...
  'Units', 'normalized', ...
  'HandleVisibility','callback', ...
  'Position', [0.25 0.87 0.7 0.11]);
box(hPlotAxes1)

hPlotAxes2 = axes('Parent', hFigure, ...
  'Units', 'normalized', ...
  'HandleVisibility','callback', ...
  'Position',[0.25 0.72 0.7 0.11]);
box(hPlotAxes2)

hPlotAxes3 = axes('Parent', hFigure, ...
  'Units', 'normalized', ...
  'HandleVisibility','callback', ...
  'Position',[0.25 0.55 0.7 0.11]);
box(hPlotAxes3)

hPlotAxes4 = axes('Parent', hFigure, ...
  'Units', 'normalized', ...
  'HandleVisibility','callback', ...
  'Position',[0.25 0.38 0.7 0.11]);
box(hPlotAxes4)

hPlotAxes5 = axes('Parent', hFigure, ...
  'Units', 'normalized', ...
  'HandleVisibility','callback', ...
  'Position',[0.25 0.21 0.7 0.11]);
box(hPlotAxes5)

hPlotAxes6 = axes('Parent', hFigure, ...
  'Units', 'normalized', ...
  'HandleVisibility','callback', ...
  'Position',[0.25 0.05 0.7 0.10]);
box(hPlotAxes6)

% edit buttons + labels for top 3 plots
h.freqEdit1 = uicontrol('Style', 'edit',...
  'units','normalized', ...
  'Position', [0.11 0.955 0.08 0.03],...
  'Callback', @freq1);

h.freqText1 = uicontrol('Style','text',...
  'units','normalized', ...
  'Position',[-0.01 0.95 0.12 0.03],...
  'FontSize', 8, ...
  'String', 'Frequency');

h.ampEdit1 = uicontrol('Style', 'edit',...
  'units','normalized', ...
  'Position', [0.11 0.91 0.08 0.03],...
  'Callback', @amp1);

h.ampText1 = uicontrol('Style','text',...
  'units','normalized', ...
  'Position',[-0.01 0.905 0.12 0.03],...
  'FontSize', 8, ...
  'String','Amplitude');

h.phaseEdit1 = uicontrol('Style', 'edit',...
  'units','normalized', ...
  'Position', [0.11 0.865 0.08 0.03],...
  'Callback', @phase1);

h.phaseText1 = uicontrol('Style','text',...
  'units','normalized', ...
  'Position',[-0.01 0.86 0.12 0.03],...
  'FontSize', 8, ...
  'String','Phase');

h.freqEdit2 = uicontrol('Style', 'edit',...
  'units','normalized', ...
  'Position', [0.11 0.805 0.08 0.03],...
  'Callback', @freq2);

h.ampEdit2 = uicontrol('Style', 'edit',...
  'units','normalized', ...
  'Position', [0.11 0.76 0.08 0.03],...
  'Callback', @amp2);

h.phaseEdit2 = uicontrol('Style', 'edit',...
  'units','normalized', ...
  'Position', [0.11 0.715 0.08 0.03],...
  'Callback', @phase2);

h.freqEdit3 = uicontrol('Style', 'edit',...
  'units','normalized', ...
  'Position', [0.11 0.635 0.08 0.03],...
  'Callback', @freq3);

h.ampEdit3 = uicontrol('Style', 'edit',...
  'units','normalized', ...
  'Position', [0.11 0.59 0.08 0.03],...
  'Callback', @amp3);

h.phaseEdit3 = uicontrol('Style', 'edit',...
  'units','normalized', ...
  'Position', [0.11 0.545 0.08 0.03],...
  'Callback', @phase3);

% noise slider
h.noise4 = uicontrol('Style', 'slider',...
  'units','normalized', ...
  'Min', 0, 'Max',10, ...
  'Value', 0, ...
  'SliderStep', [0.1,0.1]/10,...
  'Position', [0.06 0.4 0.10 0.02],...
  'Callback', @noise);

h.noise4text = uicontrol('Style','text',...
  'FontUnits', 'normalized',...
  'units','normalized', ...
  'Position',[0.05 0.45 0.1 0.02],...
  'String', 'Noise');

% filter
h.filter = uicontrol('Style', 'edit',...
  'units','normalized', ...
  'Position', [0.06 0.24 0.08 0.03],...
  'Callback', @filter);

h.filterText = uicontrol('Style','text',...
  'units','normalized', ...
  'FontUnits', 'normalized',...
  'Position', [0.025 0.29 0.15 0.02],...
  'String','LP Filter');

% plot power/phase
h.outputSelection = uibuttongroup('Visible', 'on',...
  'units','normalized', ...
  'Position', [0.05 0.06 0.125 0.075],...
  'SelectionChangeFcn', @outputSelection);

h.outputSelection1 = uicontrol(h.outputSelection,  'Style', 'radiobutton',...
  'units','normalized', ...
  'FontUnits', 'normalized',...
  'Position', [0.05 0.5 1 0.4],...
  'FontSize', 1, ...
  'String','Power');

h.outputSelection2 = uicontrol(h.outputSelection,  'Style', 'radiobutton',...
  'units','normalized', ...
  'FontUnits', 'normalized',...
  'Position', [0.05 0 1 0.4],...
  'FontSize', 1, ...
  'String','Phase');

% set defaults
h.sigDur      = 3;
h.sampRate    = 2000;
h.sampPeriod  = 1/h.sampRate;
h.time        = 0:h.sampPeriod:h.sigDur-h.sampPeriod;

h.freq1  =  0;
h.amp1   =  1;
h.phase1 =  0;

h.freq2  =  0;
h.amp2   =  1;
h.phase2 =  0;

h.freq3  =  0;
h.amp3   =  1;
h.phase3 =  0;

h.noise = 0;
sig4 = zeros(size(h.time));

lowPassFilter = false;

power = true;
phase = false;

plotSignal

% turn on GUI
set(hFigure,'Visible','on')
set(hFigure, 'HandleVisibility', 'on')

% signal frequency
  function freq1(hObject, ~)
    h.freq1 = str2num(get(hObject, 'String'));
    if isempty(h.freq1)
      h.freq1 = 0;
    end
    plotSignal
  end

  function freq2(hObject, ~)
    h.freq2 = str2num(get(hObject, 'String'));
    if isempty(h.freq2)
      h.freq2 = 0;
    end
    plotSignal
  end

  function freq3(hObject, ~)
    h.freq3 = str2num(get(hObject, 'String'));
    if isempty(h.freq3)
      h.freq3 = 0;
    end
    plotSignal
  end

  function amp1(hObject, ~)
    h.amp1 = str2num(get(hObject, 'String'));
    if isempty(h.amp1)
      h.amp1 = 0;
    end
    plotSignal
  end

  function amp2(hObject, ~)
    h.amp2 = str2num(get(hObject, 'String'));
    if isempty(h.amp2)
      h.amp2 = 0;
    end
    plotSignal
  end

  function amp3(hObject, ~)
    h.amp3 = str2num(get(hObject, 'String'));
    if isempty(h.amp3)
      h.amp3 = 0;
    end
    plotSignal
  end

  function phase1(hObject, ~)
    h.phase1 = str2num(get(hObject, 'String'));
    if isempty(h.phase1)
      h.phase1 = 0;
    end
    plotSignal
  end

  function phase2(hObject, ~)
    h.phase2 = str2num(get(hObject, 'String'));
    if isempty(h.phase2)
      h.phase2 = 0;
    end
    plotSignal
  end

  function phase3(hObject, ~)
    h.phase3 = str2num(get(hObject, 'String'));
    if isempty(h.phase3)
      h.phase3 = 0;
    end
    plotSignal
  end

  function noise(~, ~)
    h.noise = get(h.noise4,'value');
    plotSignal
  end

  function filter(hObject, ~)
    h.filterValue = str2num(get(hObject, 'String'));
    if isempty(h.filterValue)
      lowPassFilter = false;
    else
      lowPassFilter = true;
    end
    plotSignal
  end

  function outputSelection(~, ~)
    power = get(h.outputSelection1, 'value');
    phase = get(h.outputSelection2, 'value');
    plotSignal
  end


  function plotSignal
    
    maxFreq = max([h.freq1, h.freq2, h.freq3]);
    maxAmp  = max([h.amp1, h.amp2, h.amp3 max(abs(sig4))]);
    
    sig1 = sin(2*h.freq1*h.time*pi + h.phase1) * h.amp1;
    cla(hPlotAxes1)
    plot(hPlotAxes1, h.time, sig1, 'b');
    if maxAmp ~= 0
      set(hPlotAxes1, 'YLim', [-maxAmp maxAmp]);
    end
    
    sig2 = sin(2*h.freq2*h.time*pi + h.phase2) * h.amp2;
    cla(hPlotAxes2)
    plot(hPlotAxes2, h.time, sig2, 'b');
    if maxAmp ~= 0
      set(hPlotAxes2, 'YLim', [-maxAmp maxAmp]);
    end
    
    sig3 = sin(2*h.freq3*h.time*pi + h.phase3) * h.amp3;
    cla(hPlotAxes3)
    plot(hPlotAxes3, h.time, sig3, 'b');
    if maxAmp ~= 0
      set(hPlotAxes3, 'YLim', [-maxAmp maxAmp]);
    end
    
    if h.noise == 0
      sig4 = zeros(size(sig1));
    else
      sig4 = 0 + h.noise.*randn(size(sig1));
    end
    cla(hPlotAxes4)
    plot(hPlotAxes4, h.time, sig4, 'b');
    if maxAmp ~= 0
      set(hPlotAxes4, 'YLim', [-maxAmp maxAmp]);
    end
    
    % calculate combined signal
    sig5 = sig1+sig2+sig3+sig4;
    cla(hPlotAxes5)
    plot(hPlotAxes5, h.time, sig5, 'b');
    if max(sig5) ~= 0
      set(hPlotAxes5, 'YLim', [-max(abs(sig5)) max(abs(sig5))]);
    end
    
    % axis labels
    xlabel('Time (S)', 'Parent', hPlotAxes5)
    ylabel('Amplitde', 'Parent', hPlotAxes5)
    
    if lowPassFilter
      hold(hPlotAxes5,'on')
      sig6 = ft_preproc_lowpassfilter(sig5, h.sampRate, h.filterValue, 6, 'but', 'twopass', false);
      plot(hPlotAxes5, h.time, sig6, 'r', 'LineWidth', 2);
    else
      sig6 = [];
    end
    
    % plot power
    cla(hPlotAxes6)
    f1  = fft(sig5)/length(h.time);
    hz = linspace(0, h.sampRate/2, length(h.time)/2+1);
    
    % power or phase
    xlabel('Frequency (Hz)', 'Parent', hPlotAxes6)
    if power
      stem(hPlotAxes6, hz, (abs(f1(1:length(f1)/2+1))*2).^2)
      ylabel('Power', 'Parent', hPlotAxes6)
    elseif phase
      stem(hPlotAxes6, hz, angle(f1(1:length(f1)/2+1)))
      ylabel('Phase', 'Parent', hPlotAxes6)
      set(hPlotAxes6, 'ytick',-pi:pi/2:pi)
    end
    
    if lowPassFilter
      f2  = fft(sig6)/length(h.time);
      hold(hPlotAxes6,'on')
      if power
        stem(hPlotAxes6, hz, (abs(f2(1:length(f2)/2+1))*2).^2, 'r')
        ylabel('Power', 'Parent', hPlotAxes6)
      elseif phase
        stem(hPlotAxes6, hz, angle(f2(1:length(f2)/2+1)), 'r')
        ylabel('Phase', 'Parent', hPlotAxes6)
        set(hPlotAxes6, 'ytick',-pi:pi/2:pi)
      end
    end
    
    if maxFreq ~= 0
      set(gca,'xlim',[0 maxFreq+0.2*maxFreq])
    end
    
  end

end
