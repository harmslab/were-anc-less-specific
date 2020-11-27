__description__ = \
"""
Helper code for converting raw values from plate reader into fluorescence
anisotropy binding experiments.
"""
__author__ = "Michael J. Harms"
__date__ = "2020-09-01"

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import scipy.optimize

import re

def calculate_r(vv,vh,G=1.0):
    """
    Calculate anisotropy from vertical (vv), horizontal (vh), and empirical
    G-factor (G). 
    """

    return (vv - G*vh)/(vv + 2*G*vh)

def species(Mt,Xt,Kd):
    """
    Return species of M, X, and MX for a single-site binding model given
    total protein concentration (Mt), total ligand concentration (Xt),
    and dissociation constant (Kd).
    """

    a = 1
    b = -(Xt + Mt + Kd)
    c = Xt*Mt

    roots = [(-b + np.sqrt(b**2 - 4*a*c))/(2*a),
             (-b - np.sqrt(b**2 - 4*a*c))/(2*a)]

    allowed_roots = []
    for r in roots:
        if np.isreal(r):
            if r >= 0 and r <= Mt and r <= Xt:
                allowed_roots.append(r)

    if len(allowed_roots) == 0:
        err = "no root found!\n"
        raise RuntimeError(err)
    if len(allowed_roots) > 1:
        err = "no unique root found!\n"
        raise RuntimeError(err)

    MX = allowed_roots[0]
    X = Xt - MX
    M = Mt - MX

    return M,X,MX

def get_peptide_Kd_scalar(Kd,Mt,Xt):
    """
    returns factor to multiply Kd_app by to obtain actual Kd.

    Kd: probe Kd
    Mt: total protein conc
    Xt: total probe conc

    Follows Eq. 7 in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3000649/
    """

    M, X, MX = species(Mt,Xt,Kd)
    delta_0 = MX/Xt
    peptide_Kd_scalar = 1 + Xt*(1-delta_0/2)/Kd+delta_0/(1-delta_0)

    return peptide_Kd_scalar

def read_plate_layout(plate_layout):
    """
    Read a plate layout.  See the read_file doc string for details of the
    format.

    Returns a list, where each element is a plate, that has the form:

    [{"protein":{(row_name,column_number):protein,
      "peptide":{(row_name,column_number):peptide,
      "value":{(row_name,column_number):value},
      ...,]

    If the contents of a cell can be coerced into a float, value will be a
    float.  If the contents of a cell is a string that starts with na (case-
    insensitive), the cell will be np.nan.  Otherwise, value will be a string.
    """

    plates = [{"protein":{},"peptide":{},"conc":{}}]
    plate_counter = 0
    df = pd.read_excel(plate_layout,header=None)
    for i in range(len(df)):
        row = np.array(df.iloc[i])
        if row[0] == "protein":
            column_protein = row[:]
            continue
        elif row[0] == "peptide":
            column_peptide = row[:]
            continue
        else:
            try:
                np.isnan(row[0])
                plate_counter += 1
                plates.append({"protein":{},"peptide":{},"conc":{}})
                continue
            except TypeError:
                row_name = row[0]

        for j in range(1,len(row)):
            plates[-1]["protein"][(row_name,j)] = column_protein[j]
            plates[-1]["peptide"][(row_name,j)] = column_peptide[j]
            try:
                plates[-1]["conc"][(row_name,j)] = np.float(row[j])
            except ValueError:
                plates[-1]["conc"][(row_name,j)] = row[j]

    return plates


def read_file(plate_file,plate_layout,min_value=5e4,G=1.0):
    """
    plate_file: output from the plate reader
    plate_layout: excel file denoting plate layout
    min_value: if not none, make any value below this NaN.  only applies to
               values in the main plate, not any named data.
    G: G-factor for anisotropy

    Inputs
    ------

    Expects plate_layout to be a pandas readable spreadsheet with the form that
    follows.  For a plate with "M" rows and "N" columns:

    protein,prot1,prot2,prot3,...,protN
    peptide,pep1,pep2,pep3,...,pepN
    A,concA1,concA2,concA3,...,concAN
    B,concB1,concB2,concB3,...,concBN
    ...
    M,concM1,concM2,concM3,...,concMN

    where prot1 is the protein in column 1, pep3 is the peptide in column 3,
    etc.  This assumes a column only has one protein/peptide pair present.

    It is also possible to have named data.  If you were to, for example, put
    the string "blank" in a cell, this will be read out as named_data.  (see
    below).

    If the contents of a cell (in plate_layout) can be coerced into a float, the
    cell will be read as a concentration.  If it is a string that starts with
    na (case-insensitive) it will be read as a concentration of value np.nan
    (probably an empty cell).  Otherwise, the cell will be interpreted as a key
    for some named data.

    Returns
    -------

    A data frame and a dictionary.

    Data frame holds:
    plate, plate_row, plate_column, c1, c2, protein, peptide, and concentration.

    Dictionary holds any named cells in a plate. For example, if plate 0 had two
    cells named blank with c1 and c2 of [10,11] and [44,55], this would return:

    {
        0:{
            "blank":{
                "c1":np.array([10,11]),
                "c2":np.array([44,55]),
                }
            }
    }

    """

    # Read the plate layout file
    plates = read_plate_layout(plate_layout)

    # Read plate_file dumped by plate reader
    f = open(plate_file,'rb')
    data = f.read().decode("utf16")
    f.close()

    # Create out dictionary for main data
    out_dict = {"plate":[],
                "plate_row":[],
                "plate_column":[],
                "c1":[],
                "c2":[]}

    # Dictionary to hold named data
    named_data = {}

    num_plates = 0
    plate_counter = 0
    next_row = "header"
    for l in data.split("\r\n"):
        columns = l.split("\t")

        # Get number of plates if no plates have been read
        if num_plates == 0:
            if len(columns) > 0 and columns[0].startswith("##BLOCKS"):
                num_plates = int(columns[0].split("=")[1])
                continue

        # Read a plate
        if num_plates > 0 and len(columns) > 1 and columns[0] == "":

            # real_columns holds columns with data
            real_columns = columns[2:-1]

            # Parse header
            if real_columns[0] == "A1":

                 # Sanity check to make sure this row looks like it should
                if next_row != "header":
                    err = f"expecting {next_row} row, but got header row.\n"
                    raise ValueError(err)

                plate_counter += 1
                out_dict["plate"].extend([int(plate_counter)
                                          for _ in range(len(real_columns))])
                out_dict["plate_row"].extend([r[0] for r in real_columns])
                out_dict["plate_column"].extend([int(r[1:]) for r in real_columns])
                next_row = "c1"

            else:

                # Sanity check to make sure this row looks like it should
                if next_row not in ["c1","c2"]:
                    err = f"expecting {next_row} row, but got data row.\n"
                    raise ValueError(err)

                # Parse columns
                for c in real_columns:
                    try:
                        value = float(c)
                    except ValueError:
                        value = np.nan

                    out_dict[next_row].append(value)

                # Update next row counter
                if next_row == "c1":
                    next_row = "c2"
                else:
                    next_row = "header"

    # Sanity check for plates
    if num_plates != len(plates):
        err = "number of plates different in data file and plate layout.\n"
        raise ValueError(err)

    # Make data frame from output
    df = pd.DataFrame(out_dict)

    # Now use plate_layout to map these values to main data frame or
    # named data.
    conc = []
    protein = []
    peptide = []
    for i in range(len(df)):

        # Get plate, row, and column
        p = df.iloc[i].plate
        row = df.iloc[i].plate_row
        col = df.iloc[i].plate_column
        c = plates[p-1]["conc"][(row,col)]

        # If this is a concentration, record it
        if type(c) is float:
            conc.append(c)
            protein.append(plates[p-1]["protein"][(row,col)])
            peptide.append(plates[p-1]["peptide"][(row,col)])

        # If this is not a concentration, it's going to be named data
        else:

            # make conc, prot, peptide nan
            conc.append(np.nan)
            protein.append(np.nan)
            peptide.append(np.nan)

            # Have we seen named data for this plate?
            try:
                named_data[p]
            except KeyError:
                named_data[p] = {}

            # Have we seen this name before on this plate?
            try:
                named_data[p][c]
            except KeyError:
                named_data[p][c] = {}
                named_data[p][c]["c1"] = []
                named_data[p][c]["c2"] = []

            # Record the two channels
            named_data[p][c]["c1"].append(df.iloc[i].c1)
            named_data[p][c]["c2"].append(df.iloc[i].c2)

    # Update dictionary
    df["protein"] = protein
    df["peptide"] = peptide
    df["conc"] = conc

    # Set things below min value to nan
    if min_value is not None:
        mask = np.logical_or(df.c1 < min_value,df.c2 < min_value)
        df.loc[mask,"c1"] = np.nan
        df.loc[mask,"c2"] = np.nan

    # Remove nan values
    df = df[np.logical_not(np.isnan(df.conc))]
    df = df[np.logical_not(np.isnan(df.c1))]
    df = df[np.logical_not(np.isnan(df.c2))]
    df["r"] = calculate_r(df.c2,df.c1,G=G)

    # Clean up
    df["plate"] = df["plate"].astype(int)
    df["plate_column"] = df["plate_column"].astype(int)

    # Convert named_data into numpy arrays
    for p in named_data.keys():
        for c in named_data[p].keys():
            named_data[p][c]["c1"] = np.array(named_data[p][c]["c1"])
            named_data[p][c]["c2"] = np.array(named_data[p][c]["c2"])
            named_data[p][c]["r"] = calculate_r(named_data[p][c]["c2"],
                                                named_data[p][c]["c1"],
                                                G=G)
    return df, named_data

def average_tech_reps(df,
                      remove_outlier=True,
                      outlier_min_std=0.01,
                      min_next_ratio=0.1,
                      same_dist_ratio=0.7):
    """
    Take the dataframe output from read_file and average the tehnical 
    replicates.  Do some data clean up.

    df: data frame from read_file
    If remove_outlier is True, the function will look for tech reps with a 
    standard deviation greater than "outlier_min_std."  An outlier will be 
    identified if the ratio between the distances between the furthest points
    is less than min_next_ratio and the ratio between the distances between 
    the closest two points is grater than same_dist_ratio. 
  
    outlier
         +                               +
          +
    
    no outlier

        +                 +              +

    no outlier

        +++
    """

    # Prep for output data frame
    fit_dict = {"protein":[],
                "peptide":[],
                "plate":[],
                "plate_row":[],
                "plate_column":[],
                "conc":[],
                "r":[],
                "r_err":[],
                "weight":[],
                "outlier_removed":[]}

    # mask of obs to keep
    to_keep = [np.array((0,1),dtype=np.int),
               np.array((0,2),dtype=np.int),
               np.array((1,2),dtype=np.int)]

    # Go through all concentrations
    for c in np.unique(df.conc):

        # Get values for this concentration 
        df_c = df[df.conc == c]
        columns = np.array(df_c.plate_column)

        # See if there is an outlier
        outlier_removed = False
        v = np.array(df_c.r)
        if  len(v) == 3 and np.std(v) > outlier_min_std:

            D = np.array(((v[0] - v[1])**2,(v[0] - v[2])**2,(v[1] - v[2])**2))
            closest = to_keep[np.argmin(D)]

            D.sort()
            if D[0]/D[1] < min_next_ratio and D[1]/D[2] > same_dist_ratio:
                v = v[closest]
                columns = columns[closest]
                outlier_removed = True

        # record to write out in data frame
        fit_dict["protein"].append(df_c["protein"].iloc[0])
        fit_dict["peptide"].append(df_c["peptide"].iloc[0])
        fit_dict["conc"].append(c*1e-6)
        fit_dict["r_err"].append(np.std(v))
        fit_dict["r"].append(np.mean(v))
        fit_dict["weight"].append(1.0) 
        fit_dict["outlier_removed"].append(outlier_removed)
        fit_dict["plate"].append(df_c["plate"].iloc[0])
        fit_dict["plate_row"].append(df_c["plate_row"].iloc[0])
        fit_dict["plate_column"].append(columns)

    fit_df = pd.DataFrame(fit_dict)

    return fit_df

def load_data(experiments,
              remove_outlier=True,
              peptides=["A5cons",
                        "A6cons",
                        "phage_ctl_0",
                        "phage_ctl_1",
                        "phage_ctl_2",
                        "phage_ctl_4",
                        "phage_ctl_5",
                        "phage_ctl_6",
                        "phage_ctl_7",
                        "phage_ctl_8",
                        "phage_ctl_9"]):
    """
    Convenience function that allows one to load a whole bunch of experiments,
    with different peptides, into a single data frame.

    experiments should be a list of dictionaries of the following form:

            [{"protein":"hA6",
              "name_in_file":"hA6_4.3",
              "Kd":45,
              "prot_conc":4.2,
              "probe_conc":4.2,
              "data_file":"13_main-collection.txt",
              "plate_file":"13_plate-layout.xlsx"},...]

    remove_outlier: whether or not to look for outlier points and remove them
                    when averaging technical reps
    peptides: list of peptides.  these are used to build regular expressions
              to match peptides in each data file.  It looks for an exact match
              at the start of the string, allowing any trailing characters.
              NOTE: this could lead to problems if you had peptides with names
              like pep10, pep100.
    """


    pep_patterns = [re.compile(f"{p}") for p in peptides]

    proteins = set([e["protein"] for e in experiments])

    times_pep_was_seen = dict([(protein,dict([(p,0) for p in peptides]))
                               for protein in proteins])

    all_df = []
    for expt in experiments:
        df, _ = read_file(expt["data_file"],expt["plate_file"])
        df = df[df.protein == expt["name_in_file"]]

        peptide_Kd_scalar = get_peptide_Kd_scalar(Kd=expt["Kd"],
                                                  Mt=expt["prot_conc"],
                                                  Xt=expt["probe_conc"])

        peps_in_df = np.unique(df.peptide)
        for p in peps_in_df:

            for pattern in pep_patterns:
                if pattern.match(p):

                    pep_df = df[df.peptide == p]
                    plates = np.unique(pep_df.plate)

                    protein = expt["protein"]
                    peptide = pattern.pattern
                    for plate in plates:

                        times_pep_was_seen[protein][peptide] += 1

                        single_rep = pep_df[pep_df.plate == plate]
                        fit_df = average_tech_reps(single_rep,remove_outlier=remove_outlier)

                        fit_df["protein"] = protein
                        fit_df["peptide"] = peptide
                        fit_df["rep_number"] = times_pep_was_seen[protein][peptide]
                        fit_df["Kd_scalar"] = peptide_Kd_scalar

                        fit_df["plate_file"] = expt["plate_file"]
                        fit_df["data_file"] = expt["data_file"]
                        fit_df["name_in_file"] = expt["name_in_file"]
                        fit_df["plate_number"] = plate

                        all_df.append(fit_df)

                    break


    return pd.concat(all_df)

# -----------------------------------------------------------------------------
# Functions for fitting binding models
# -----------------------------------------------------------------------------

def single_site(x,Kd,A,B):
    """
    First-order binding model model. 
    """
    
    return A + B*(x/(x + Kd))

def single_site_r(param,x,obs,weights,fixed_A,fixed_B):
    """
    Residuals function for first-order model.
    """
    
    Kd = param[0]
    if fixed_A is None:
        A = param[1]
    else:
        A = fixed_A
        
    if fixed_B is None:
        B = param[2]
    else:
        B = fixed_B
    
    return (single_site(x,Kd,A,B) - obs)*weights

def fit_model(x,obs,weights,param_guesses=(1,1,1),fixed_A=None,fixed_B=None):
    """
    Fit the first-order model.
    """

    if ((fixed_A is None or fixed_B is None) and not (fixed_A is None and fixed_B is None)):
        err = "You must fix neither A nor B, or both A and B.\n"
        raise ValueError(err)
    
    fit = scipy.optimize.least_squares(single_site_r,
                                       param_guesses,
                                       args=(x,obs,weights,fixed_A,fixed_B))
    fit_Kd = fit.x[0]
    if fixed_A is not None:
        fit_A = fixed_A
    else:
        fit_A = fit.x[1]
    if fixed_A is not None:
        fit_B = fixed_B
    else:
        fit_B = fit.x[2]
    
    return fit_Kd, fit_A, fit_B


def fit_and_plot(df_list,
                 baseline_list=None,
                 offset_to_reference=False,
                 required_change_in_signal=None,
                 name_list=None,
                 color_list=None,
                 plot=True,
                 title=None,Kd_guess=10e-6,
                 xlim=None,ylim=None,
                 fit_fontsize=16,
                 keep_fit=True,log=False,
                 xlabel="[peptide] (uM)",ylabel="anisotropy",
                 alpha=1,legend=True,plot_err=True,
                 fig=None,ax=None):

    if plot and fig is None:
        fig, ax = plt.subplots(1,1,figsize=(8,8))

    all_fits = []
    for i, fit_df in enumerate(df_list):

        A_guess = np.max(fit_df.r)
        B_guess = np.min(fit_df.r) - A_guess
        
        if baseline_list is not None:
            fixed_A = baseline_list[i][0]
            fixed_B = baseline_list[i][1]
        else:
            fixed_A = None
            fixed_B = None
        
        try:
            fit_Kd, fit_A, fit_B = fit_model(fit_df.conc,
                                             fit_df.r,
                                             fit_df.weight,
                                             param_guesses=(Kd_guess,A_guess,B_guess),
                                             fixed_A=fixed_A,
                                             fixed_B=fixed_B)
            conditions = [fit_Kd > 1e-7,
                          fit_Kd < 1e-3]
            
            if required_change_in_signal is not None:
                conditions.append(fit_B < required_change_in_signal)

            if np.sum(conditions) != len(conditions):
                keep_fit = False
                
        except ValueError:
            keep_fit = False
            
        if not keep_fit:
            fit_Kd = np.nan
            fit_A = np.nan
            fit_B = np.nan
            
        all_fits.append((fit_Kd,fit_A,fit_B))
        
        # Offset data so it starts at 0
        if keep_fit and offset_to_reference:
            fit_df.loc[:,"r"] = fit_df.loc[:,"r"] - fit_A
            fit_Kd, fit_A, fit_B = fit_model(fit_df.conc,
                                             fit_df.r,
                                             fit_df.weight,
                                             param_guesses=(Kd_guess,A_guess,B_guess),
                                             fixed_A=0,fixed_B=fit_B)
            
        if not plot:
            return None, None, all_fits
            
        
        kwargs = {"fmt":"o",
                  "capsize":4,
                  "alpha":alpha}

        if color_list is not None:
            kwargs["color"] = color_list[i]
        if name_list is not None:

            if np.isnan(fit_Kd):
                label = "{} (Kd: NaN uM)".format(name_list[i])
            else:
                label = "{} (Kd: {:.1f} uM)".format(name_list[i],fit_Kd*1e6)

            kwargs["label"] = label
            

        if plot_err:
            err = fit_df.r_err
        else:
            err = None
            
        ax.errorbar(fit_df.conc*1e6,fit_df.r,err,**kwargs)

        if keep_fit:

            if log:
                min_x = np.min(fit_df.conc)/10
                max_x = np.max(fit_df.conc)*10
                xrange = np.exp(np.linspace(np.log(min_x),np.log(max_x)))
            else:
                min_x = 0
                max_x = np.max(fit_df.conc)*1.1
                xrange = np.linspace(min_x,max_x,100)

            kwargs = {"lw":3,"alpha":alpha}
            if color_list is not None:
                kwargs["color"] = color_list[i]

            ax.plot(xrange*1e6,single_site(xrange,fit_Kd,fit_A,fit_B),**kwargs)

    if xlim is not None:
        ax.set_xlim(*xlim)

    if ylim is not None:
        ax.set_ylim(*ylim)

    if log:
        ax.set_xscale("log")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if title is not None:
        ax.set_title(title)

    if name_list is not None and legend:
        ax.legend()

    return fig, ax, all_fits
