
from urn_packages import plt

#Function to add labels to plot
def add_labels(plt,myxlabel,myylabel,font_sizes,leg):
    plt.xlabel(myxlabel,fontsize=font_sizes[0])
    plt.ylabel(myylabel,fontsize=font_sizes[1])
    plt.grid(color='grey', linestyle='--', linewidth=0.5)
    plt.xticks(fontsize=font_sizes[2])
    plt.yticks(fontsize=font_sizes[3])
    if(leg==1):
        plt.legend(fontsize=font_sizes[4])
    return plt
