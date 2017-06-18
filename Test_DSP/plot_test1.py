
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

#% Demo macro plot 4 bars and give a different color to each one.
#% Also plots the value of the bar above the bar.
clc
#% Clear the command window.
plt.close(all)
#% Close all figures (except those of imtool.)
clear
#% Erase all existing variables.
workspace
#% Make sure the workspace panel is showing.
fontSize = 30.
format(compact)
#% Ask user for the number of bars.
defaultValue = 4.
titleBar = 'Enter an integer value'
userPrompt = 'Enter the number of bars'
caUserInput = inputdlg(userPrompt, titleBar, 1., cellarray(np.hstack((num2str(defaultValue)))))
if isempty(caUserInput):
    return []


#% Bail out if they clicked Cancel.
integerValue = np.round(str2double(cell2mat(caUserInput)))
#% Check for a valid integer.
if np.isnan(integerValue):
    #% They didn't enter a number.  
#% They clicked Cancel, or entered a character, symbols, or something else not allowed.
integerValue = defaultValue
message = sprintf('I said it had to be an integer.\nI will use %d and continue.', integerValue)
uiwait(warndlg(message))

#% Define sample data in the range 20-80.
x = np.arange(1., (integerValue)+1)
y = 20.+80.*np.random.rand(integerValue)
numberOfBars = length(y)
button = menu('Use which colormap?', 'Custom', 'Random', 'Jet', 'Hot', 'Lines')
if button == 1.:
    #% Make up a custom colormap specifying the color for each bar series.
barColorMap[0,:] = np.array(np.hstack((.2, .71, .3)))
#% Green Color for segment 1.
barColorMap[1,:] = np.array(np.hstack((.25, .55, .79)))
#% Blue Color for segment 2.
barColorMap[2,:] = np.array(np.hstack((.9, .1, .14)))
#% Red Color for segment 3.
barColorMap[3,:] = np.array(np.hstack((.9, .9, .14)))
#% Yellow Color for segment 4.
#% I have not defined any more than 4 colors in this demo.
#% For any number of bars beyond 4, just make up random colors.
if numberOfBars > 4.:
    barColorMap[4:numberOfBars,0:3.] = np.random.rand((numberOfBars-4.), 3.)


elif button == 2.:
    #% Example of using colormap with random colors
    barColorMap = np.random.rand(numberOfBars, 3.)
    
elif button == 3.:
    #% Example of using pre-defined jet colormap
    barColorMap = plt.jet(numberOfBars)
    
elif button == 4.:
    #% Example of using pre-defined Hot colormap
    barColorMap = plt.hot(numberOfBars)
    
else:
    #% Example of using pre-defined lines colormap
    barColorMap = lines(numberOfBars)
    

#% Plot each number one at a time, calling bar() for each y value.
for b in np.arange(1., (numberOfBars)+1):
    #% Plot one single bar as a separate bar series.
    
#% Fancy up the graph.
plt.grid(on)
caption = sprintf('Data plotted in %d barseries, each with a different color', length(y))
plt.title(caption, 'FontSize', fontSize)
plt.xlabel('x', 'FontSize', fontSize)
plt.ylabel('y', 'FontSize', fontSize)
#% Restore the x tick marks.
set(plt.gca, 'XTickMode', 'Auto')
#% set(gca, 'XTickLabels', xTickLabels);
#% Enlarge figure to full screen.
set(plt.gcf, 'units', 'normalized', 'outerposition', np.array(np.hstack((0., 0., 1., 1.))))
#% Give a name to the title bar.
set(plt.gcf, 'name', 'Demo by ImageAnalyst', 'numbertitle', 'off')