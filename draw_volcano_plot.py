# Modified from kevincoakley metabolomicsworkbench/jupyter-notebooks/MWPerformVolcanoPlotAnalysis.ipynb

import numpy as np
import seaborn as sns

def DrawVolcanoPlot(VolcanoPlotDataFrame, condition1, condition2, LogFoldChangeColID = "log2(FoldChange)", PValueColID = "P-value", LogPValueColID ="-log10(P-value)", LogFoldChangeThreshold = 1.0, PValueThreshold = 0.05, minlogPValue = 0.05, PlotStyle = "darkgrid", FontScale = 1.3, Title = "Volcano Plot", TitleFontWeight = "bold", LabelsFontWeight = "bold", PlotWidth = 9, PlotHeight = 6):
    
    # Setup color for data points...
    DataFrame = VolcanoPlotDataFrame
    
    # Significant up foldchange...
    ColorColID = 'Color'
    DataFrame.loc[(DataFrame[LogFoldChangeColID] >= LogFoldChangeThreshold) 
                  & (DataFrame[PValueColID] <= PValueThreshold), ColorColID] = "r"
    # Significant down foldchange...
    DataFrame.loc[(DataFrame[LogFoldChangeColID] <= -LogFoldChangeThreshold) 
                  & (DataFrame[PValueColID] <= PValueThreshold), ColorColID] = "b"
    # Significant change...
    DataFrame.loc[(DataFrame[LogFoldChangeColID] >= -LogFoldChangeThreshold)
                  & (DataFrame[LogFoldChangeColID] <= LogFoldChangeThreshold)
                  & (DataFrame[PValueColID] <= PValueThreshold), ColorColID] = "grey"
    # Intermediate change...
    DataFrame[ColorColID].fillna('grey', inplace = True)
    
    ColorsPalette = {"r" : "r", "b" : "b", "grey" : "grey"}
    
    sns.set(rc = {'figure.figsize':(PlotWidth, PlotHeight)})
    sns.set(style = PlotStyle, font_scale = FontScale)
    
    Axis = sns.scatterplot(x = LogFoldChangeColID, y = LogPValueColID, hue = ColorColID, data = VolcanoPlotDataFrame,
                           palette = ColorsPalette)

    # Draw vertical lines at LogFoldChangeThreshold...
    Axis.axvline(-LogFoldChangeThreshold, color = 'violet', linewidth = 1, linestyle = "dashed")
    Axis.axvline(LogFoldChangeThreshold, color = 'violet', linewidth = 1, linestyle = "dashed")
    
    # Draw a horizontal line at -log10(PValueThreshold)
    HLinePos = -np.log10(PValueThreshold)
    Axis.axhline(HLinePos, color = 'green', linewidth = 1, linestyle = "dashed")
    
    # Set title and labels...
    Axis.set_title(Title, fontweight = TitleFontWeight)
    Axis.set_xlabel(LogFoldChangeColID, fontweight = LabelsFontWeight)
    Axis.set_ylabel(LogPValueColID , fontweight = LabelsFontWeight)
    
    # Set figure legent...
    handles, labels  =  Axis.get_legend_handles_labels()
    Axis.legend(handles, ['Not significant', 'Significantly increased in '+ condition2, 'Significantly decreased in '+ condition2], bbox_to_anchor=(1.04,1), loc="upper left")
    
    # top metabolites in significant up foldchange
    top_significant_foldchange = DataFrame.loc[(DataFrame[LogFoldChangeColID] >= LogFoldChangeThreshold) 
                  & (DataFrame[LogPValueColID] >= -np.log10(minlogPValue))]
    
    top_significant_foldchange = top_significant_foldchange.sort_values(LogPValueColID, ascending = False)
    
    return top_significant_foldchange


def DrawVolcanoPlotColored(VolcanoPlotDataFrame, condition1, condition2, class_level_col, color_palette, LogFoldChangeColID = "log2(FoldChange)", PValueColID = "P-value", LogPValueColID ="-log10(P-value)", LogFoldChangeThreshold = 1.0, PValueThreshold = 0.05, minlogPValue = 0.05, PlotStyle = "darkgrid", FontScale = 1.3, Title = "Volcano Plot", TitleFontWeight = "bold", LabelsFontWeight = "bold", PlotWidth = 9, PlotHeight = 6):
      
    # Setup color for data points by class level indicator...
    DataFrame = VolcanoPlotDataFrame.replace('nan', np.nan).dropna(subset=[class_level_col])
    unique_classes = DataFrame[class_level_col].unique()
    class_color = dict(zip(unique_classes,color_palette))
    class_color['Insignificant'] = 'grey'

    DataFrame['color'] = DataFrame[class_level_col].map(class_color)

    palette = {i:i for i in DataFrame.color.unique()}
    palette['grey']='grey'
    
    # Significant up foldchange...
    ColorColID = 'color'

    # Significant change...
    DataFrame.loc[(DataFrame[LogFoldChangeColID] >= -LogFoldChangeThreshold)
                  & (DataFrame[LogFoldChangeColID] <= LogFoldChangeThreshold)
                  & (DataFrame[PValueColID] <= PValueThreshold), ColorColID] = "grey"

    # Intermediate change...
    DataFrame.loc[(DataFrame[PValueColID] >= PValueThreshold), ColorColID] = "grey"

    sns.set(rc = {'figure.figsize':(PlotWidth, PlotHeight)})
    sns.set(style = PlotStyle, font_scale = FontScale)

    Axis = sns.scatterplot(x = LogFoldChangeColID, y = LogPValueColID, hue = ColorColID, data = DataFrame,
                           palette = palette)

    # Draw vertical lines at LogFoldChangeThreshold...
    Axis.axvline(-LogFoldChangeThreshold, color = 'violet', linewidth = 1, linestyle = "dashed")
    Axis.axvline(LogFoldChangeThreshold, color = 'violet', linewidth = 1, linestyle = "dashed")

    # Draw a horizontal line at -log10(PValueThreshold)
    HLinePos = -np.log10(PValueThreshold)
    Axis.axhline(HLinePos, color = 'green', linewidth = 1, linestyle = "dashed")

    # Set title and labels...
    Axis.set_title(Title, fontweight = TitleFontWeight,y=1.10)
    Axis.set_xlabel(LogFoldChangeColID, fontweight = LabelsFontWeight)
    Axis.set_ylabel(LogPValueColID , fontweight = LabelsFontWeight)

    # Set figure legent...
    handles, labels = Axis.get_legend_handles_labels()

    newLabels, newHandles = [], []
    for handle, label in zip(handles, labels):
        for key_class_name in class_color:
            if class_color[key_class_name] == label:
                newLabels.append(key_class_name)
                newHandles.append(handle)

    Axis.legend(newHandles,newLabels,bbox_to_anchor=(1.04,1), loc="upper left")

    xmin, xmax, ymin, ymax = Axis.axis()

    Axis.annotate('', xy=(LogFoldChangeThreshold, ymax), xycoords='data', xytext=(xmax-(xmax/4), ymax), 
                arrowprops=dict(arrowstyle="<-", color='r', lw=4))

    Axis.text(LogFoldChangeThreshold, ymax+(ymax/25), 'Increased in '+condition2, fontsize=14)

    Axis.annotate('', xy=(-LogFoldChangeThreshold, ymax), xycoords='data', xytext=(xmin-(xmin/4), ymax), 
                arrowprops=dict(arrowstyle="<-", color='b', lw=4))

    Axis.text(-LogFoldChangeThreshold, ymax+(ymax/25), 'Decreased in '+condition2, fontsize=14, ha="right")
    
    return Axis
