#!/usr/bin/env python

################################################################################
# multi-analysis.py

# purpose: This code plots the results of the code multi.
# Author: Nicolas Baillot d'Etivaux

################################################################################

# Imported packages:
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec 
import matplotlib.ticker as ticker
from matplotlib import cm
from matplotlib.offsetbox import AnchoredText
import numpy as np
import netCDF4 as nc

# Allow the use of tex format for the labels:
plt.rcParams.update({"text.usetex": True, "font.size" : 15})

################################################################################
def netcdf_extract(res_dir):
    """ Extract the netcdf data from the result directory res_dir:
        Args: 
            res_dir: result directory in which is stored the netcdf output file.
        Return: 
            ds: netCDF4 dataset containing the results.
    """
    ds=nc.Dataset(res_dir+"/output.nc","r")
    return ds
################################################################################
################################################################################
def lanczos_vs_planczosif_plot(ds,res_dir,outer_iterations):
    """ Plot the J, Jb and Jo:
        Args: 
            res_dir (string): result directory in which is stored the netcdf output file.
        Return:
            (void): Plots the cost functions and save them as out_name_J,Jb or Jo.
    """
    lanczos,planczosif={},{}
    keys=['j_nl','jo_nl','jb_nl','j','jo','jb','rho_sqrt','beta']
    labels=[r'$J^{nl}$',r'$J_o^{nl}$',r'$J_b^{nl}$',r'$J$',r'$J_o$',r'$J_b$',r'$\sqrt{\rho}$',r'$\beta$']
    for k in keys:
        lanczos[k]=[]
        planczosif[k]=[]
    
    for outer in ds.groups:
        for k in keys:
            lanczos[k].append(np.array(ds[outer]['lanczos'][k][:]))
            planczosif[k].append(np.array(ds[outer]['planczosif'][k][:]))
            
    diff={}    
    for ik,k in enumerate(lanczos.keys()):
        diff[k]=[]
        lanczos[k]=np.reshape(lanczos[k],-1)
        planczosif[k]=np.reshape(planczosif[k],-1)
        for i in range(len(lanczos[k])):
            diff[k].append(lanczos[k][i]-planczosif[k][i])
    
        # Plot the results:
        out_name=res_dir+'/compare_{}.png'.format(k)
        print("plotting:", out_name)
        ylabel12=labels[ik]
        ylabel3='diff'
        xmax=max(len(lanczos[k]),len(planczosif[k]))
        x=range(xmax)
        xlabel='iterations'
        legend=['Lanczos','PlanczosIF']
        compare_plots(lanczos[k],planczosif[k],ylabel12,diff[k],ylabel3,x,xlabel,outer_iterations,legend,out_name)
################################################################################
################################################################################
def obs_plot(ds,res_dir):
    """Plots the observations:
    """
    out_name=res_dir+'/obs.png'
    x_obs=np.array(ds['x_obs'][:])
    y_obs=np.array(ds['y_obs'][:])
    obs_val=np.array(ds['obs_val'][:])
    
    print("plotting:", out_name)
    fig = plt.figure()
    plt.scatter(x_obs,y_obs,c=obs_val,cmap=plt.get_cmap('copper'))
    plt.colorbar()
    plt.savefig(out_name)
    plt.clf()
    plt.close()
################################################################################
################################################################################
def coord_plot(ds,res_dir):
    """Check the model grid (should plot a straight line x=y):
    """
    for io in ds.groups:
        out_name=res_dir+'/grid_coord_{}.png'.format(io)
        x_coord=np.array(ds[io]['x_coord'][:])
        y_coord=np.array(ds[io]['y_coord'][:])
        fig = plt.figure()
        plt.scatter(x_coord,y_coord)
        plt.savefig(out_name)
        plt.clf()
        plt.close()
################################################################################
################################################################################
def hxg_plot(ds,res_dir):
    """Plots hxg
    """
    x_obs=np.array(ds['x_obs'][:])
    y_obs=np.array(ds['y_obs'][:])
    for io in ds.groups:
        for algo in ds[io].groups:
            hxg=np.array(ds[io][algo]['hxg'][:])
            out_name=res_dir+'/'+algo+'_hxg_'+io+'.png'
            print("plotting:", out_name)
            fig=plt.figure()
            plt.scatter(x_obs,y_obs,c=hxg,cmap=plt.get_cmap('copper'))
            plt.colorbar()
            plt.savefig(out_name)
            plt.clf()
            plt.close()
################################################################################
################################################################################
def innovation_plot(ds,res_dir):
    """Plots the innovation
    """
    x_obs=np.array(ds['x_obs'][:])
    y_obs=np.array(ds['y_obs'][:])
    for io in ds.groups:
        for algo in ds[io].groups:
            d=np.array(ds[io][algo]['d'][:])
            out_name=res_dir+'/'+algo+'_d_'+io+'.png'
            print("plotting:", out_name)
            fig=plt.figure()
            plt.scatter(x_obs,y_obs,c=d,cmap=plt.get_cmap('copper'))
            plt.colorbar()
            plt.savefig(out_name)
            plt.clf()
            plt.close()
################################################################################
################################################################################
def compare_plots(obj1,obj2,ylabel12,obj3,ylabel3,x,xlabel,outer_iterations,legend,out_name):
    """Produces usefull comparision plots bteween obj1 and obj2:
    Args:
        out_name (string): Name of the output png file.
        obj1, obj2 (list): The two lists of numbers to compare.
        yabel12 (string) : Label of the two objects obj1 and obj2.
        obj3 (list)      : List of numbers regarding which one wants to compare
                           obj1 and obj2 (for example obj1-obj2).
        ylabel3 (string) : Label of obj3.
        x (list)         : List of numbers for x-axis.
        xlabel (string)  : Label of x-axis.
        legend (list(string): Labels for obj1 and obj2.
    Return:
        (void): Creates the plot and save it in the png file named by out_name arg.
    """
    print("plotting:",out_name)
    # Create figure window to plot data
    fig = plt.figure(1, figsize=(9,9))
    gs = gridspec.GridSpec(2, 1, height_ratios=[6, 3])

    # Top plot: obj1 vs obj2
    ax1 = fig.add_subplot(gs[0])
    ax1.plot(x[:len(obj1)],obj1,color='black',label=legend[0])
    ax1.plot(x[:len(obj2)],obj2,color='black',linestyle='dashed',label=legend[1])
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
               ncol=2, mode="expand", borderaxespad=0.)

    plt.subplots_adjust(hspace=0.5)
    ymax=max(max(obj1),max(obj2))
    if ymax > 100:
        plt.yscale("log")
    ax1.set_ylabel(ylabel12)
    start, end = ax1.get_xlim()
    ax1.vlines(outer_iterations, 0, ymax, colors='blue', linestyles='dashed')

    # Bottom plot: obj3
    ax2 = fig.add_subplot(gs[1])
    ax2.plot(x[:len(obj3)],obj3,color='black')
    ax2.axhline(color="gray", zorder=-1)
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel(ylabel3)
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    plt.savefig(out_name)
    plt.clf()
    plt.close()
################################################################################
################################################################################
def field_plot(matrix,out_name):
    """Represents matrix using matshow.
    """
    print('plotting:',out_name)
    fig=plt.figure()
    plt.matshow(np.array(matrix),cmap=plt.get_cmap('copper'))
    plt.colorbar()
    plt.savefig(out_name)
    plt.clf()
    plt.close()
################################################################################
################################################################################
def field_plot2(matrix,x_coord,y_coord,out_name):
    """Represents matrix using scatterplot.
    """
    print('plotting:',out_name)
    fig=plt.figure()
    val=[]
    for i in range(len(x_coord)):
        val.append(matrix[i][j])        
    plt.scatter(x_coord,y_coord,c=val,cmap=plt.get_cmap('copper'))
    plt.colorbar()
    plt.savefig(out_name)
    plt.clf()
    plt.close()
################################################################################
################################################################################
def compare_methods_plot(compare_methods_data,methods_list,outer_iterations,compare_methods_dir):
    """Plots the comparision between the different methods:
    """
    lanczos,planczosif,obj_list,diff_list={},{},{},{}
    keys=['j_nl','jo_nl','jb_nl','j','jo','jb','rho_sqrt','beta']
    labels=[r'$J^{nl}$',r'$J_o^{nl}$',r'$J_b^{nl}$',r'$J$',r'$J_o$',r'$J_b$',r'$\sqrt{\rho}$',r'$\beta$']
    legend=[]
    for met in methods_list:
        legend_met=[]
        for space in ['control','model']:
            legend_met.append(met+'-'+space)
        legend.append(legend_met)
    
    for k,key in enumerate(keys):
        lanczos[key]=[]
        planczosif[key]=[]
        obj_list[key]=[]
        diff_list[key]=[]
        for i,ds in enumerate(compare_methods_data):
            lanczos[key].append([])
            planczosif[key].append([])
            for io in ds.groups:
                lanczos[key][i].append(np.array(ds[io]['lanczos'][key][:]))
                planczosif[key][i].append(np.array(ds[io]['planczosif'][key][:]))
            lanczos[key][i]=np.reshape(lanczos[key][i],-1)
            planczosif[key][i]=np.reshape(planczosif[key][i],-1)
            obj_list[key].append([lanczos[key][i],planczosif[key][i]])
            # Compute the difference between lanczos and planczosif:
            diff_tmp=[]
            for j in range(len(lanczos[key][i])):
                diff_tmp.append(lanczos[key][i][j]-planczosif[key][i][j])
            diff_list[key].append(diff_tmp)
        ylabel1=labels[k]
        ylabel2="diff"
        xmax=0
        for obj in obj_list[key]:
            xmax=max(len(obj[0]),len(obj[1]))
        x=list(range(xmax))
        xlabel="iterations"
        out_name=compare_methods_dir+'/compare_{}.png'.format(key)
        print('plotting:',out_name)
        compare_plots_2N(obj_list[key],ylabel1,diff_list[key],ylabel2,x,xlabel,outer_iterations,legend,out_name)
################################################################################
################################################################################
def compare_methods_2D_outer(compare_methods_data,methods_list,compare_methods_dir):
    """Plots the comparision between the different methods:
    """
    diff_dict={}
    keys=['xg']
    labels=[r'$x^g$']
    # legend=[]
    # for met in methods_list:
    #     legend_met=[]
    #     for space in ['control','model']:
    #         legend_met.append(met+'-'+space)
    #     legend.append(legend_met)

    ds_th=compare_methods_data[0]
    for io in ds_th.groups:
        for algo in ds_th[io].groups:
            diff_dict[algo]={}
            for key in keys:
                diff_dict[algo][key]=[]
                for m,met in enumerate(methods_list):
                    ds=compare_methods_data[m]
                    data_to_compare=np.array(ds[io][algo][key][:])
                    diff_dict[algo][key].append(np.array(data_to_compare))
                # Compute the difference as a 2D matrix -- diff = (method - theoretical)
                for m,met in enumerate(methods_list[1:]):
                    diff_matrix=np.zeros(np.shape(diff_dict[algo][key][m]))
                    for line in range(len(diff_dict[algo][key][m])):
                        for col in range(len(diff_dict[algo][key][m][line])):
                            diff=diff_dict[algo][key][m][line][col]-diff_dict[algo][key][0][line][col]
                            diff_matrix[line][col]=diff
                    # Plot the difference:
                    out_name=compare_methods_dir+'/theoretical_vs_{}/{}_{}_{}.png'.format(met,algo,key,io)
                    print('plotting:',out_name)
                    field_plot(diff_matrix,out_name)
################################################################################
################################################################################
# /!\ Obsolete (tmp):
def vec_plot(results_file,column_of_interest,label,out_name):
    """Plots 1D vectors:
    """
    print("plotting:",out_name)
    # Get the data:
    outer_vectors=np.genfromtxt(results_file, comments='#')
    io=1
    vec=[]
    indices=[]
    vec_io=[]
    indices_io=[]
    for i  in range(len(outer_vectors)):
        # Add the elemnts of the outer vector for iteration io:
        if int(outer_vectors[i,0])==io:
            vec_io.append(outer_vectors[i,column_of_interest])
            indices_io.append(outer_vectors[i,1])
        # Add the outer vector for outer iteration io to vec:    
        else:
            io=io+1
            vec.append(vec_io)
            indices.append(indices_io)
            vec_io=[]
            indices_io=[]
            vec_io.append(outer_vectors[i,column_of_interest])
            indices_io.append(outer_vectors[i,1])
    # Add the last outer vector:
    vec.append(vec_io)
    indices.append(indices_io)
    
    if not len(vec)==1:
        fig, subplots = plt.subplots(len(vec),1)
        for i, ax in enumerate(subplots):
            plt.ticklabel_format(axis="y", style="sci", scilimits=(-3,3))
            ax.plot(indices[i][:],vec[i][:],color='blue')
            ax.set_ylabel(label+'_{}$'.format(i))
        plt.tight_layout()
        fig.align_ylabels(subplots[:])    
        plt.subplots_adjust(hspace=0.5)
    else:
        fig=plt.figure()
        plt.plot(indices[0][:],vec[0][:],color='blue')
        plt.ylabel(label+r'_{}$'.format(0))
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    plt.savefig(out_name)
    plt.clf()
    plt.close()
################################################################################
################################################################################
def compare_plots_2N(obj_list,ylabel1,diff_list,ylabel2,x,xlabel,outer_iterations,legend,out_name):
    """Produces usefull comparision with 2N plots:
    Args:
        out_name (string): Name of the output png file.
        obj_list (list)      : List containing lists to be compared together.
        yabel (string)   : Label of the y-axis.
        diff_list (list)      : List of numbers regarding which one wants to compare.
        ylabel (string)  : Label of diff_list.
        x (list)         : List of numbers for x-axis.
        xlabel (string)  : Label of x-axis
    Return:
        (void): Creates the plot and save it in the png file named by out_name arg.
    """
    color_map=cm.get_cmap('copper', len(obj_list))
    colors=color_map(range(len(obj_list)))
    # Create figure window to plot data:
    fig = plt.figure(1, figsize=(9,9))
    gs = gridspec.GridSpec(2, 1, height_ratios=[6, 2])

    # Top plot: [obj_list1,obj_list2,..obj_listN]:
    ax1 = fig.add_subplot(gs[0])
    plt.yscale("log")
    ymax=0.
    for o,obj in enumerate(obj_list):
        ymax_tmp=max(max(obj[0]),max(obj[1]))
        if (ymax_tmp>ymax):
            ymax=ymax_tmp
        ax1.plot(x[:len(obj[0])],obj[0],color=colors[o],label=legend[o][0])
        ax1.plot(x[:len(obj[1])],obj[1],color=colors[o],linestyle='dashed',label=legend[o][1])
        plt.legend(bbox_to_anchor=(0., 1.102 , 1., .102), loc='lower left',
                   ncol=3, mode="expand", borderaxespad=0.1)
        plt.subplots_adjust(top=0.8)
    ax1.set_ylabel(ylabel1)
    #start, end = ax1.get_xlim()
    #ax1.xaxis.set_ticks(np.arange(start, end, 1.))
    ax1.vlines(outer_iterations, 0, ymax, colors='blue', linestyles='dashed')

    # Bottom plot: diff_list
    ax2 = fig.add_subplot(gs[1])
    for d,diff in enumerate(diff_list):
        ax2.plot(x[:len(diff)],diff,color=colors[d])
    ax2.axhline(color="gray", zorder=-1)
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel(ylabel2)
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    #ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.4f'))
    plt.savefig(out_name)
    plt.clf()
    plt.close()
################################################################################
################################################################################
# /!\ Obsolete: (tmp)
def lmp_compare(out_names,lmp_to_compare,column_of_interest,ylabel1,ylabel2,outer_iterations_list):
    """Compare spectral and ritz lmp modes:
    """
    try:
        legend=[]
        for lmp_mode in ['ritz','spectral','none']:
            legend.append([lmp_mode+'-model',lmp_mode+'-control'])
            
        for r,res_dirs in enumerate(lmp_to_compare):
            diff_list=[]
            obj_list=[]
            for res_dir in res_dirs:
                res1=np.genfromtxt(res_dir+'PlanczosIF_model_space.dat', comments='#')    
                res2=np.genfromtxt(res_dir+'lanczos_control_space.dat', comments='#')
                diff=[]
                l_iter=min(len(res1[0]),len(res2[0]))
                for c in range(len(res1)):
                    diff_col=[]
                    for l in range(l_iter):
                        diff_col.append(res1[c,l]-res2[c,l]) 
                    diff.append(diff_col)
                diff=np.array(diff)
                diff_list.append(diff[:,column_of_interest])
                obj_list.append([res1[:,column_of_interest],res2[:,column_of_interest]])

            # maybe dirtyish but...
            itot=list(range(len(diff)))
            print("plotting lmp comparision for:\n",out_names[r],"\n for column:", column_of_interest)
            x=itot
            xlabel='iterations'
            out_name=out_names[r]
            outer_iterations=outer_iterations_list[r]
            compare_plots_2N(out_name,obj_list,ylabel1,diff_list,ylabel2,x,xlabel,outer_iterations,legend)
    except:
        print("Error with lmp comparision of:\n",res_dirs,"\n")
##################################################################################
