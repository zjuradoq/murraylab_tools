# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import string
import math
import csv
import collections
import os
import warnings
from ..utils import *
from datetime import date as pydate

__all__ = ["dna2nM_convert", "echo_round", "AbstractMixture", "WellReaction",
           "MixtureMaterial", "TXTLMasterMix", "SourcePlate", "EchoSourceMaterial",
           "Pick", "EchoRun", "DestinationPlate", "dead_volume", "max_volume",
           "usable_volume"]

dead_volume = 15000 + 6000 # Dead volume per well in an Echo source plate, in nL
max_volume  = 65000 # Maximum Echoable volume in an Echo source plate, in nL
usable_volume = max_volume - dead_volume # Maximum pipettable volume per well
                                         # in an Echo source plate, in nL

def dna2nM_convert(dnaconc, dnalength):
    '''
    Converts DNA ng/ul to nM

    DNA_conc(ng/ul) * 1e6 (uL/L) * (bp*mol)/660g * 1/dna_length(bp = DNA (nM)
    Double stranded DNA is 660g/(bp*mol), single stranded DNA is 330g/(bp*mol)
    '''
    return (dnaconc*1e6)/(660*dnalength)

def echo_round(x, base = 25):
    '''
    Round a volume to an exact number the Echo can pipette

    Arguments:
        x - Desired volume.
        base - Size of an Echo drop, in dimensionality of the volume (default nL).

    Returns: x rounded to the nearest multiple of base.
    '''
    return int(base * round(float(x)/base))

def floatify(element):
    '''
    Convert a string to a float, if possible; otherwise, return the string
    right back. Code ripped off code from Jacob Gabrielson's answer on
    StackOverflow thread 736043.

    Empty strings are converted to "0" to make vector addition work on relevant
    ranges.
    '''
    if not type(element) == str:
        raise ValueError("Why are you trying to floatify a " + type(element) + \
                         "?")
    if element == "":
        return 0.0
    partition = element.partition('.')
    if element.isdigit() or \
       (partition[0].isdigit() and partition[1]=='.' and \
        partition[2].isdigit()) \
        or \
        (partition[0]=='' and partition[1]=='.' and partition[2].isdigit()) \
        or \
        (partition[0].isdigit() and partition[1]=='.' and partition[2]==''):
        return float(element)
    else:
        return element

def process_column_argument(col):
    '''
    Convenience function for processing a column-name argument. Converts
    from either a string ("C") or a zero-indexed int (2) to a zero-indexed
    int. Nonetype arguments are returned as None.
    '''
    if col == None:
        return None
    if type(col) == str:
        if col == "":
            raise ValueError("Column argument can't be an empty string.")
        upper_col = col.upper()
        col_num = 0
        while len(upper_col) > 0:
            col_num *= 26
            col_num += string.ascii_uppercase.find(upper_col[-1]) + 1
            upper_col = upper_col[:-1]
        #Remember, it's zero-indexed!
        col_num -= 1
        return col_num
    elif type(col) == int:
        return col
    else:
        raise TypeError("Column argument must be a string ('C') or a " +\
                        "zero-indexed int (2)")


class Pick():
    '''
    Container class for making a single Echo pick from some source well to
    a destination well.
    '''
    def __init__(self, source_material, source_well, destination_well, volume):
        self.source_material   = source_material
        self.source_well       = source_well
        self.destination_well  = destination_well
        self.volume            = int(volume)


class EchoSourceMaterial():
    '''
    Container class for a single material on a source plate. Keeps track of how
    much material has been used and where they need to go. When requested,
    will commit to wells (assigned by a SourcePlate) and send a list of picks.
    '''
    def __init__(self, name, concentration, length = 0, units = None):
        '''
        name -- A string to associate with the name of this material.
        concentration -- the concentration of this material, in units given by
                            "units".
        length -- Length of the DNA in bp, if material is double-stranded DNA.
                  Otherwise, length should be 0. If a (nonzero) length is given,
                  then units are assumed to be ng/uL, and are immediately
                  converted to nM.
        units -- Name of the units of this material, e.g. "nM", "x", or "ng/uL".
                    This argument is ignored if length != 0. This is only for
                    human readability purposes, and does not affect
                    calculations.

        '''
        self.name   = name
        self.length = length
        self.plate  = None
        if units is None:
            self.units = ""
        else:
            self.units  = units

        # Check for commas in the name.
        if "," in self.name:
            warnings.warn(("Material %s has comma in its name; this may bug " +
                          "the echo when you run it.") % self.name, Warning)

        # DNA concentration in ng/uL, or other reagent concentration in nM.
        concentration = float(concentration)
        self.concentration = concentration
        if length == None or length <= 0:
            self.nM = self.concentration
        else:
            self.nM = dna2nM_convert(concentration, length)
            self.ng_per_ul = concentration

        self.wells    = None
        self.picklist = []
        self.pipettelist = []
        self.total_volume_requested = 0
        self.echo_volume_requested  = 0
        self.well_volumes = None
        self.current_well = -1

    def __str__(self):
        if self.length > 0:
            conc_string = "ng/µL"
            conc = self.ng_per_ul
        else:
            conc_string = self.units
            conc = self.nM
        return "%s (%.3f %s)" % (self.name, conc, conc_string)

    def request_material(self, destination_well, volume, pipette_by_hand=False):
        '''
        Request some material be made available for picking.

        volume is the volume of material requested, in nL.

        Updates this material's total requested volume and list of
        destinations.
        '''
        if pipette_by_hand:
            self.pipettelist.append(Pick(self, None, destination_well, volume))
            self.total_volume_requested += volume
        else:
            actual_volume = echo_round(volume)
            if actual_volume == 0:
                warnings.warn(f"Requesting 0 volume from material {self.name}"
                              f" into well {destination_well}; are you sure "
                              "you want to do this?")
            else:
                self.total_volume_requested += actual_volume
                self.echo_volume_requested += actual_volume
                self.picklist.append(Pick(self, None, destination_well,
                                          actual_volume))

    def request_picklist(self):
        '''
        Commit to wells, and return a list of finalized picks from this
        material.
        '''
        usable_volume  = max_volume - dead_volume
        n_source_wells = math.ceil(float(self.echo_volume_requested) \
                                         / usable_volume)
        if n_source_wells == 0:
            print(("Warning: Material %s is requesting 0 wells in its " +\
                  "source plate to give %f total volume") % \
                  (self.name, self.echo_volume_requested))
            return
        if self.wells == None:
            self.wells = self.plate.request_wells(int(n_source_wells))
        if len(self.wells) < 1:
            warnings.warn(("Material %s has requested no wells. Are you sure "+\
                          "this is correct?") % self.name)

        self.well_volumes    = np.zeros(len(self.wells))
        self.well_volumes[0] = dead_volume
        self.current_well    = 0
        self.total_volume_requested += dead_volume
        self.echo_volume_requested += dead_volume

        for pick in self.picklist:
            volume_requested = pick.volume # For error-writing purposes.
            available_volume = max_volume - self.well_volumes[self.current_well]

            while pick.volume > available_volume:
                yield Pick(self, self.wells[self.current_well],
                           pick.destination_well, available_volume)
                self.well_volumes[self.current_well] = max_volume
                self.current_well += 1
                self.total_volume_requested += dead_volume
                self.echo_volume_requested += dead_volume
                if self.current_well >= len(self.wells):
                    raise ValueError(("Material %s has been asked to donate " +\
                                      "too much material with a call for " +\
                                      "%d nL") % \
                                    (self.name, volume_requested))
                self.well_volumes[self.current_well] = dead_volume
                pick.volume -= available_volume
                available_volume = usable_volume

            self.well_volumes[self.current_well] += pick.volume
            yield pick

    def text_recipe(self):
        ret_str = "\n%s:\n\tstock concentration: %.2f %s" % \
                    (self.name, self.nM, self.units)
        if self.length > 0:
            ret_str += " (%.2f ng/uL)" % self.concentration
        ret_str += "\n\ttotal volume: %.2f uL" % \
                    (self.total_volume_requested / 1000.0)
        return ret_str

class EchoRun():
    '''
    Defines and prints an Echo picklist from one of several algorithms.

    Parameters:
        -- rxn_vol: Volume of a single TX-TL reaction, in nL. Default 10000.
        -- DPtype: Destination plate type. Should be a string recognizable by
                    the Echo Plate Reformat software, though that is not
                    enforced in this code. Default "Nunc_384_black_glassbottom"
        --plate: Source plate. Should be a SourcePlate object. Required for
                    TX-TL setup experiments, but not for association spreadsheet
                    experiments.
    '''
    def __init__(self, rxn_vol = 10000, DPtype = None, plate = None,plate_log_folder=None):
        # The user gives the reaction volume in microliters; store it in nL
        self.rxn_vol = rxn_vol

        if DPtype == None:
            self.DPtype = "Nunc_384_black_glassbottom"
        else:
            self.DPtype = DPtype
        if '384' in self.DPtype:
            self.rows = 16
            self.cols = 24
        elif '96' in self.DPtype:
            self.rows = 8
            self.cols = 12
        else:
            raise ValueError(("Unrecognized plate type %s. If this is a " + \
                              "valid plate type, contact a murraylab_tools " +\
                              "developer to have it added.") % self.DPtype)
        if plate == None:
            self.plates = [SourcePlate(plate_log_folder, "384PP_AQ_BP")]
        else:
            self.plates = [plate]

        if plate_log_folder is None:
            self.plate_log_folder = ""
        else:
            self.plate_log_folder = plate_log_folder
        self.material_dict   = dict()
        self.reactions       = dict()
        self.picklist        = []
        self.fill_material   = None

    def define_plate(self, SPname, SPtype, DPtype):
        '''
        Convenience function for setting plate parameters.
        '''
        self.SPname = SPname
        self.SPtype = SPtype
        self.DPtype = DPtype

    def add_material(self, material):
        '''
        Add a material to the materials list for this EchoRun, if an identical
        material is not already in the list. Attempting to add a material with
        the same name but different properties as another material already in
        this object's material list will raise a ValueError.

        Returns 0 if the material was added successfully; returns 1 if an
        identical material was already in this EchoRun object's material list,
        so the new material was not added.
        '''
        if material.name in self.material_dict.keys():
            prior_mat = self.material_dict[material.name]
            if prior_mat.name == material.name \
               and prior_mat.concentration == material.concentration \
               and prior_mat.length == material.length \
               and prior_mat.units  == material.units \
               and prior_mat.plate == material.plate:
                return 1
            else:
                raise ValueError(
                    (f"Tried to add material {material.name} with concentration"
                     f" {material.concentration} {material.units}, length "
                     f"{material.length}, plate {material.plate}; that material"
                     f" already exists with concentration "
                     f"{prior_mat.concentration} {prior_mat.units}, length "
                     f"{prior_mat.length}, plate {prior_mat.plate}."))
        elif isinstance(material, MixtureMaterial) \
             and material.rxn_vol != self.rxn_vol:
            warnings.warn(f"EchoRun object has rxn_vol = {self.rxn_vol}, "
                             f"but was given a MixtureMaterial with rxn_vol = "
                             f"{material.rxn_vol}. Assuming EchoRun's rxn_vol"
                             " is correct.")
            material.rxn_vol = self.rxn_vol

        else:
            self.material_dict[material.name] = material
            material.plate = self.plates[0]
            return 0


    def build_picklist_from_txtl_setup_excel(self, input_filename):
        '''
        CURRENTLY NONFUNCTIONAL DO NOT USE

        Build an Echo picklist based on a TX-TL setup spreadsheet (v2.1 or
        newer).
        '''
        raise NotImplementedError("'build_picklist_from_txtl_setup_excel' " + \
                                  "hasn't been implemented yet. Use " + \
                                  "'build_picklist_from_txtl_setup_csvs'.")
        # Open the workbook and identify all the important sheets
        # workbook = pyxl.load_workbook(input_filename)
        # recipe_sheet = None
        # stock_sheet = None
        # layout_sheet = None
        # for sheet in workbook:
        #     if sheet.title == "Recipe":
        #         recipe = sheet
        #     elif sheet.title == "Stocks":
        #         stocks = sheet
        #     elif sheet.title == "Layout":
        #         layout = sheet

    def build_picklist_from_txtl_setup_csvs(self, stock_filename,
                                            recipe_filename):
        '''
        Build an Echo picklist based on a pair of CSV documents produced from
        a TX-TL setup spreadsheet (v>=2.1).

        The stock sheet is a CSV describing materials on the source plate
        (plasmids, usually).

        The recipe sheet is a CSV describing the master mix and what materials
        from the stock sheet go in what destination wells, in what quantity.

        The layout sheet is not used by this function.
        '''
        self.make_master_mix = True

        ####################
        # Read Input Files #
        ####################

        # Read in stock file
        stock_sheet = np.empty(shape = (12,5), dtype = object)
        with mt_open(stock_filename, 'rU') as stock_file:
            stock_reader = csv.reader(stock_file)
            rownum = -1
            for row in stock_reader:
                rownum += 1
                for colnum in range(len(row)):
                    element = floatify(row[colnum])
                    stock_sheet[rownum, colnum] = element

        # Read in recipe file
        recipe_sheet = np.zeros(shape = (384+20, 16), dtype = object)
        with mt_open(recipe_filename, 'rU') as recipe_file:
            recipe_reader = csv.reader(recipe_file)
            rownum = -1
            for row in recipe_reader:
                rownum += 1
                if rownum >= recipe_sheet.shape[0]:
                    print("Warning -- You are trying to add to more than " + \
                          "384 wells in the destination plate. " + \
                          "Extra wells will be clipped.")
                    break
                for colnum in range(len(row)):
                    element = floatify(row[colnum])
                    if element:
                        recipe_sheet[rownum, colnum] = element

        # Set some magic numbers based on the recipe file
        self.rxn_vol = float(recipe_sheet[11,10]) * 1e3
        self.extract_fraction = float(recipe_sheet[12,2])
        self.buffer_fraction  = 0.75 - self.extract_fraction
        self.mm_excess = float(recipe_sheet[11,12])


        ######################
        # Register Materials #
        ######################

        # Magic numbers here are spreadsheet positions of volumes to add.
        material_total_vols  = [0]*10
        for i in range(len(material_total_vols)):
            for j in range(20,recipe_sheet.shape[0]):
                if recipe_sheet[j, 5+i]:
                    material_total_vols[i] += recipe_sheet[j, 5+i] * 1e3

        # Assign source wells and register materials

        # Register Water
        self.add_material(EchoSourceMaterial("water", 1, 0, self.plates[0]))
        water = self.material_dict["water"]

        # Register other materials
        stocks = []
        for i in range(len(material_total_vols)):
            if stock_sheet[i+2,1] == "":
                for j in range(i+1, len(material_total_vols)):
                    if stock_sheet[i+2,1] != "":
                        raise RuntimeWarning("You left a blank row in your " + \
                            "stock file. This will cause alignment shifts in "+\
                            "your recipe sheet and you will almost surely " +  \
                            "use the wrong amounts of ingredients. Are you " + \
                            "sure this is what you want?")
                continue
            material_name          = stock_sheet[i+2, 1]
            if isinstance(material_name, float):
                continue
            material_concentration = stock_sheet[i+2, 2]
            material_length        = stock_sheet[i+2, 3]
            new_material           = EchoSourceMaterial(material_name,
                                                        material_concentration,
                                                        material_length,
                                                        self.plates[0])
            is_duplicate_material = self.add_material(new_material)
            if not is_duplicate_material:
                stocks.append(new_material)

        ##################
        # Register picks #
        ##################

        # Make and register the TXTL master mix.
        mm_materials = []
        for i in range(11,17):
            if recipe_sheet[i,4] is None or recipe_sheet[i,4] == 0:
                continue
            name  = recipe_sheet[i,0]
            stock = recipe_sheet[i,1]
            final = recipe_sheet[i,2]
            material = EchoSourceMaterial(name, stock, 0, self.plates[0])
            mm_materials.append((material, final))

        extract_fract = None
        buffer_fract  = None
        # This loop must be backwards, because we're removing elements as we go.
        for i in range(len(mm_materials)-1, -1, -1):
            mat_tuple = mm_materials[i]
            name = mat_tuple[0].name.lower()
            if "extract" in name:
                del mm_materials[i]
                extract_fract = mat_tuple[1]
                if extract_fract != self.extract_fraction:
                    warnings.warn(("EchoRun object and TX-TL spreadsheet have ",
                                   "different extract fraction values; ",
                                   "defaulting to spreadsheet value."))
            elif "buffer" in name:
                del mm_materials[i]
                buffer_fract = mat_tuple[1]
        if extract_fract is None:
            raise ValueError(("TX-TL spreadsheet has no recognizable extract ",
                             "in master mix."))
        if buffer_fract is None:
            raise ValueError(("TX-TL spreadsheet has no recognizable buffer ",
                             "in master mix."))

        if not "txtl_mm" in self.material_dict:
            # self.material_dict['txtl_mm'] = TXTLMasterMix(
            self.add_material(TXTLMasterMix(extract_fraction = extract_fract,
                                txtl_fraction = extract_fract + buffer_fract,
                                mm_excess = self.mm_excess,
                                rxn_vol = self.rxn_vol))
        txtl = self.material_dict["txtl_mm"]
        for material, final in mm_materials:
            txtl.add_material(material, final,
                              dimensionality = "final_concentration")

        first_row = 20
        last_row  = 20 + 384
        n_rxns    = 0
        for rownum in range(first_row, last_row):
            # Check to see if there's a well used in this row; if not, skip it.
            if recipe_sheet[rownum, 1] == 0:
                continue
            n_rxns += 1
            well = recipe_sheet[rownum, 1]
            if well == 0:
                raise ValueError(("Error on row for ID #%d of recipe sheet: " +\
                                 "Must have a destination well.") % \
                                 (rownum - 21))
            if well in self.reactions.keys():
                raise ValueError("Well %s already has a reaction!" \
                                 % well)
            self.reactions[well] = WellReaction(self.rxn_vol, well)

            # Material picks (magic number warning -- magic numbers define
            # positions of relevant blocks in the recipe sheet)
            for mat_num in range(len(material_total_vols)):
                colnum = mat_num + 5
                volume = recipe_sheet[rownum, colnum] * 1e3
                if volume != None and volume > 0:
                    source_material = stocks[mat_num]
                    self.reactions[well].add_material(source_material, volume,
                                                      dimensionality = "volume")

            # Water picks (magic number warning -- magic number defines
            # positions of relevant blocks in the recipe sheet)
            volume = recipe_sheet[rownum, 3] * 1e3
            self.reactions[well].fill_with(water)

            # Master Mix picks (magic number warning -- magic number defines
            # positions of relevant blocks in the recipe sheet)
            volume = recipe_sheet[rownum, 4] * 1e3
            self.reactions[well].add_material(txtl, volume,
                                              dimensionality = "volume")


    def load_source_plate(self, input_filename, name_col, conc_col, len_col,
                          well_col, units_col, header = True):
        '''
        Enter new materials from a CSV spreadsheet.

        Args:
            input_filename -- name of the CSV.
            name_col -- Name of the column containing the name of each material,
                        either as a string ("C") or a 0-indexed int (2)
            conc_col -- Name of the column containing the concentration of
                        each material, in ng/uL if the material is dsDNA or in
                        nM or relative concentration otherwise, either as a
                        string ("C") or a 0-indexed int (2)
            len_col  -- Name of the column containing the length of any dsDNA
                        material, either as a string ("C") or a 0-indexed int
                        (2), or None if no such column exists (for sheets
                        containing only non-dsDNA materials)
            well_col -- Name of the column containing the well location of each
                        material, either as a string ("C") or a 0-indexed int
                        (2)
            units_col -- Name of the column containing the name of the units
                        of whatever's in the well (e.g. "nM", "x", "ng/uL").
            header -- True iff there is a header row. Decides whether or not to
                        skip the first line of each file
        '''
        #####################
        # Process arguments #
        #####################
        name_idx  = process_column_argument(name_col)
        conc_idx  = process_column_argument(conc_col)
        len_idx   = process_column_argument(len_col)
        well_idx  = process_column_argument(well_col)
        units_idx = process_column_argument(units_col)

        #############
        # Read file #
        #############
        with mt_open(input_filename, 'rU') as input_file:
            reader = csv.reader(input_file)
            # Skip first row if it's a header
            if header:
                next(reader)
            for row in reader:
                name          = row[name_idx]
                concentration = floatify(row[conc_idx])
                well          = row[well_idx]
                if len_idx != None:
                    length = int(floatify(row[len_idx]))
                else:
                    length = 0
                if units_col == None:
                    units = "x"
                else:
                    units = row[units_idx]
                material = EchoSourceMaterial(name, concentration,
                                              length, units)
                self.add_material(material)
                if self.material_dict[name].wells == None:
                    self.material_dict[name].wells = [well]
                else:
                    self.material_dict[name].wells.append(well)


    def build_picklist_from_association_spreadsheet(self, input_filename,
                                                    well_column, header = True,
                                                    fill_with_water = False,
                                                    water_name = None):
        '''
        Make an Echo picklist based on an association spreadsheet, which is a
        CSV file where 1) each line is a reaction and 2) lines contains
        alternating columns of material names and final concentrations (in nM),
        starting with material name.

        Args:
            input_filename -- name of the CSV.
            well_col -- Name of the column containing the destination well of
                        each reaction, either as a string ("C") or a 0-indexed
                        int (2)
            header -- True iff there is a header row. Decides whether or not to
                        skip the first line of each file
            fill_with_water -- If true, will fill all wells to the reaction size
                                with water. Requires water_name argument.
            water_name -- Determines the name of wells containing water. Must
                            match the name given in an association spreadsheet,
                            or otherwise defined.
        '''
        well_idx = process_column_argument(well_column)
        if fill_with_water:
            if not water_name:
                raise Exception("If 'Fill with water' option selected, must "
                                "set water_name to the name of water given in "
                                "your source file.")

        with mt_open(input_filename, 'rU') as input_file:
            reader = csv.reader(input_file)
            # Skip the first row if it's a header
            if header:
                next(reader)
            for row in reader:
                if well_idx >= len(row):
                    raise ValueError("Well column out of bounds for row '%s'" %\
                                     str(row))
                well = row[well_idx]
                if well not in self.reactions:
                    self.reactions[well] = WellReaction(self.rxn_vol, well)
                i = 0
                is_name_col = True
                volume_left = self.rxn_vol
                while i < len(row):
                    # Ignore it if it's the well column.
                    if i == well_idx:
                        i += 1
                        continue
                    # First column of each pair is the material
                    if is_name_col:
                        source_material = self.material_dict[row[i]]
                        is_name_col     = False
                    # Second column of each pair is the final concentration of
                    # material
                    else:
                        final_conc = float(row[i])
                        self.reactions[well].add_material(source_material,
                                                          final_conc)
                        is_name_col = True
                    i += 1
                if fill_with_water:
                    water = self.material_dict[water_name]
                    self.reactions[well].fill_with(water)

    def dilution_series(self, material_1,       material_2,
                              material_1_final, material_2_final, first_well):
        '''Makes an NxN dilution matrix with two inputs, and adds them to a
        contiguous block on the destination plate.

        Does NOT add anything else -- in particular, no master mix or water
        are added to any wells.

        Args:
            material_1, material_2 -- First and second materials to serialy
                                        dilute (EchoSourceMaterials)
            material_1_final, material_1_final -- Lists (or numpy array)
                                                    defining the final
                                                    concentrations of each
                                                    material.
            first_well -- The upper-left-most corner of the block on the
                          destination plate you want to spit into.
        '''
        material_1.plate = self.plates[0]
        material_2.plate = self.plates[0]

        self.add_material(material_1)
        self.add_material(material_2)

        # matrix size -- check to make sure it isn't going off the plate.
        n_material_1 = len(material_1_final)
        n_material_2 = len(material_2_final)
        first_row = string.ascii_uppercase.find(first_well[0])
        first_col = int(first_well[1:])
        if first_col + n_material_2 - 1 > self.cols \
            or first_row + n_material_1 > self.rows:
            raise ValueError(("Dilution series of size %dx%d starting in " +   \
                              "well %s runs off the edge of plate of type %s") \
                             % (n_material_1, n_material_2, first_well,
                                self.DPtype))

        # Fill in matrix with picks.
        for i in range(n_material_1):
            for j in range(n_material_2):
                # Initialize the reaction.
                well = string.ascii_uppercase[first_row + i] + \
                              str(first_col + j)
                if well not in self.reactions:
                    self.reactions[well] = WellReaction(self.rxn_vol, well)
                # Diluted material picks
                self.reactions[well].add_material(material_1,
                                                  material_1_final[i])
                self.reactions[well].add_material(material_2,
                                                  material_2_final[j])

    def add_material_to_well(self, material, well, final_conc = None,
                             pipette_by_hand = False):
        '''
        Add a single material, at a single concentration, to a single well.

        Parameters:
            material - An EchoSourceMaterial object representing the material
                        to add.
            final_conc - The final concentration of material, in the same
                            units as the material
            well - Name of the well to add to.
            pipette_by_hand - If True, the material will have to
                                be added by hand by the user (and instructions
                                will be printed to that effect). Default False,
                                in which case the material will be pipetted by
                                the Echo.
        '''
        self.add_material(material)
        material.plate = self.plates[0]

        if well not in self.reactions \
           or self.reactions[well] == None:
            self.reactions[well] = WellReaction(self.rxn_vol,
                                                       well)
        if final_conc is None:
            if not isinstance(material, MixtureMaterial):
                raise ValueError(f"Tried to add non-MixtureMaterial material",
                                 f"{material} of class {type(material)} to ",
                                 f"well {well} without giving a final ",
                                 "concentration.")
            vol = material.current_vol_per_rxn()
        else:
            vol = final_conc * (self.rxn_vol / material.nM)
        self.reactions[well].add_material(material, vol,
                                          dimensionality = "volume",
                                          pipette_by_hand = pipette_by_hand)

    def add_material_to_block(self, material,  top_left,
                              bottom_right, final_conc = None,
                              pipette_by_hand = False):
        '''
        Add a single material, at a single concentration, to every well in a
        block on the destination plate.

        Parameters:
            material - An EchoSourceMaterial object representing the material
                        to add.
            final_conc - The final concentration of material, in the same units
                            as the material.
            top_left - top-left-most well of the block to add to.
            bottom_right - bottom-right-most well of the block to add to.
            pipette_by_hand - If True, the material will have to
                                be added by hand by the user (and instructions
                                will be printed to that effect). Default False,
                                in which case the material will be pipetted by
                                the Echo.
        '''
        start_row = string.ascii_uppercase.find(top_left[0])
        end_row   = string.ascii_uppercase.find(bottom_right[0])
        start_col = int(top_left[1:])-1
        end_col   = int(bottom_right[1:])-1

        for row in range(start_row, end_row+1):
            for col in range(start_col, end_col+1):
                destination = string.ascii_uppercase[row] + str(col+1)
                self.add_material_to_well(material, destination, final_conc,
                                          pipette_by_hand)

    def add_material_to_all(self, material, final_conc = None,
                            pipette_by_hand = False):
        '''
        Add a single material, at a single concentration, to every well that
        has been used so far by this EchoRun object. Use for master mixes.

        To fill all wells to maximum volume with a material, use
        'fill_all_wells_with'.

        Parameters:
            material - An EchoSourceMaterial object representing the material
                        to add.
            final_conc - The final concentration of material, in the same units
                            as the material.
            pipette_by_hand - If True, the material will have to
                                be added by hand by the user (and instructions
                                will be printed to that effect). Default False,
                                in which case the material will be pipetted by
                                the Echo.
        '''
        for well in self.reactions:
            self.add_material_to_well(material, well, final_conc,
                                      pipette_by_hand)


    def fill_well_with(self, material, well, pipette_by_hand = False):
        '''
        Adds a material to a single destination well up to the well's maximum
        volume. The material is not immediately added, just marked for later --
        you can safely add other materials to the same well after calling this
        function, and the volume of the fill material will be dynamically
        updated accordingly.

        Parameters:
            material - An EchoSourceMaterial object representing the material to
                        add.
            well - String describing which well of the destination plate to fill.
            pipette_by_hand - If True, the material will have to
                                be added by hand by the user (and instructions
                                will be printed to that effect). Default False,
                                in which case the material will be pipetted by
                                the Echo.
        '''
        self.add_material(material)
        if well not in self.reactions or self.reactions[well] == None:
            self.reactions[well] = WellReaction(self.rxn_vol, well)
        self.reactions[well].fill_with(material, pipette_by_hand)

    def fill_all_wells_with(self, material, pipette_by_hand = False):
        '''
        Same as fill_well_with, but applies to every well that has been used so
        far by this EchoRun object.
        '''
        self.add_material(material)
        for well in self.reactions:
            # self.reactions[well].fill_with(material, pipette_by_hand)
            self.fill_well_with(material, well, pipette_by_hand)

    def initialize_source_plate(self):
        """
        Initializes source plate and returns a dictionary of wells to fill:
        {(material name, material conc)--> [(well, volume)]}
        """
        material_well_dict = {}
        for mat_name in self.material_dict:
            mat = self.material_dict[mat_name]
            if mat:
                name, conc = mat.name, mat.concentration
                if mat.plate is None:
                    raise Exception(f"Material '{mat}' doesn't have a plate " +\
                                     "at plate initialization time.")
                wells_to_fill = mat.plate.request_source_wells(mat)

                material_well_dict[(name, conc)] = wells_to_fill

                if (name, conc) in mat.plate.materials_to_add:
                    material_well_dict[(name, conc)] += \
                                          mat.plate.materials_to_add[name, conc]
        #Empty Materials to_add
        materials_to_add = {}
        return material_well_dict

    def generate_picklist(self):
        for mat_name in self.material_dict:
            mat = self.material_dict[mat_name]
            if mat:
                picks = mat.plate.request_picklist(mat)
                #picks = mat.request_picklist()
                for pick in picks:
                    yield pick

    def write_picklist(self, outputname):
        '''
        Write this EchoCSVMaker's protocol to an Echo picklist, and print any
        other necessary directions for the user.

        outputname - The base filename for both the picklist and the instruction
                        file, e.g. "my_experiment/outputs/my_experiment".
        '''
        # Finalize all of the reactions.
        for reaction in self.reactions.values():
            reaction.finalize()

        # Write picklist.
        # NOTE! This MUST come before writing the comment file; comments require
        # accurate count of total_volume_requested of each material, which is
        # only calculated once all of the picks are finalized and wells are
        # committed (which happens in this block).
        with mt_open((outputname + '_EchoInput.csv'), 'w') as outcsv:
            writer = csv.writer(outcsv, lineterminator = "\n")

            # Write header
            writer.writerow(["Source Plate","Source Plate Type","Source Well",
                             "Sample ID","Destination Plate Name",
                             "Destination Well","Transfer Volume",
                             "Sample Comment"])

            #initialize source plate
            mat_well_dict = self.initialize_source_plate()
            # Write picks
            for pick in self.generate_picklist():
                if pick.source_material.name == "txtl_mm" \
                   or pick.source_material.name == "water":
                    comment = ""
                else:
                    comment = "Actual concentration: %.2f" \
                        % (pick.source_material.nM * pick.volume / self.rxn_vol)
                plate = pick.source_material.plate
                row = [plate.name, plate.type, pick.source_well,
                       pick.source_material.name, self.DPtype,
                       pick.destination_well, int(pick.volume), comment]
                writer.writerow(row)
        # Write comment file
        with mt_open((outputname+'_experiment_overview.txt'), 'w') as text_file:
            text_file.write("Materials used:\n")
            for material in self.material_dict.values():
                text_file.write(material.text_recipe())
            # Source plate loading instructions
            text_file.write("\n\nOn the source plate:")

            for material in self.material_dict.values():
                name, conc = material.name, material.concentration
                vol_list = material.well_volumes

                #Fill wells from reusable source plate
                if (name, conc) in mat_well_dict \
                        and len(mat_well_dict[name, conc]) > 0:
                    volumes_to_add = list(set(
                                     [i[1] for i in mat_well_dict[name, conc]]))
                    volumes_to_add.sort()
                    volumes_to_add.reverse()
                    wells_by_volume = \
                        {v:[i[0] for i in mat_well_dict[name, conc] \
                                 if i[1]==v] for v in volumes_to_add}
                    for vol in volumes_to_add:
                        text_file.write(f"\n\tAdd {vol/1000.0} uL of {name} in")
                        if len(wells_by_volume[vol]) > 1:
                            text_file.write(" wells: ")
                        else:
                            text_file.write(" well: ")
                        first_well = True
                        for w in wells_by_volume[vol]:
                            if first_well:
                                first_well = False
                            else:
                                text_file.write(", ")

                            text_file.write(w)
                    #for (well, vol) in mat_well_dict[name, conc]:
                    #    text_file.write("\n"+str(vol/1000.0)+"uL of "+name+" in well "+well+"\n")

            # Destination plate loading instructions (for hand-pipetted stuff).
            for material in self.material_dict.values():
                if len(material.pipettelist) > 0:
                    text_file.write("\n\nOn destination plate:")
                last_well = None

                vols = list(set([ps.volume for ps in material.pipettelist]))
                vols.sort()
                vols.reverse()
                wells_by_vol = {vol:[ps.destination_well for ps in material.pipettelist if ps.volume == vol] for vol in vols}
                for vol in wells_by_vol:

                    wells = [(w[0], int(w[1:])) for w in wells_by_vol[vol]]
                    text_file.write("\n\tPipette %.2f ul of %s into wells:" \
                                    % (vol/1000, material.name))
                    wells.sort()
                    current_row = wells[0][0]
                    for row, col in wells:
                        well = row+str(col)
                        if current_row != row:
                            text_file.write("\n\t\t"+well)
                            current_row = row
                        elif well == wells[0] and len(wells):
                            text_file.write("\n\t\t"+well)
                        else:
                            text_file.write(", "+well)

            # Make the plates write out their usage.
            for plate in self.plates:
                plate.write_to_file()

        #Write a file with the definitions (ingredients) in each well
        with mt_open((outputname + '_well_definitions.csv'), 'w') as well_defs:
            well_defs.write("well,ingredient1,concentration1,unit1,ingredient2,...\n")
            wells = list(self.reactions.keys())
            wells.sort()
            for well in wells:
                well_defs.write(well)
                for (mat, conc, unit) in self.reactions[well].materials:
                    if mat.concentration == 1:
                        unit = "fraction"
                    well_defs.write(","+mat.name+","+str(conc)+","+unit)
                well_defs.write("\n")


class AbstractMixture(object):
    '''
    Container class for mixes of liquids.

    Superclass for
        * MixtureMaterial: For generic mixes of liquids that will be used as
                            materials, including master mixes.
        * TXTLMasterMix: Subclass of MixtureMaterial specialized for TX-TL
                            master mixes.
        * WellReaction: A mix of liquids in a well on a destination plate.
    '''
    def __init__(self, vol = 0, well = None, recipe_excess = 1.0):
        self.vol       = vol
        self.finalized = False
        self.materials = [] #List of tuples (material, amount, unit)
        self.fill_material = None
        if recipe_excess < 1:
            raise ValueError("recipe_excess must be greater than or equal to 1.")
        self.recipe_excess = recipe_excess #Extra fraction (typically 1.1) of
                                           # mixture to make when printing a recipe

    def add_material(self, material, final_conc, dimensionality = "concentration"):
        '''
        Add a material at a known final concentration. Final concentrations are
        assumed to use the same dimensionality as the material (usually nM). Does NOT
        check final mix volume -- that will not be checked until
        finalize is called (happens automatically when recipe is
        called).

        Params:
            material - Usually an EchoSourceMaterial object. If using dimensionality of
                        volume or percent, material can be a string, in which
                        case a dummy EchoSourceMaterial object will be made
                        and returned for that material.
            final_conc - Final concentration of the material in this mixture.
            dimensionality - One of {"concentration", "percent", "volume"}.
                             Controls how final_conc is interpreted. If
                             dimensionality = "concentration", final_conc is the
                             final concentration in the same dimensionality as
                             the material. If "percent", will add as a fixed
                             fraction of the total reaction. If "volume", will
                             add a fixed volume, in nanoliters.
        '''
        dimensionality = dimensionality.lower()
        if dimensionality == "concentration":
            if not isinstance(material, EchoSourceMaterial):
                raise TypeError(f"Attempted to add a material with type "
                                f"{type(material)}; must be an "
                                "EchoSourceMaterial when using dimensionality of"
                                " concentration.")
            self.materials.append((material, final_conc, "concentration"))
            return

        ret_material = None
        if isinstance(material, str):
            material = EchoSourceMaterial(material, 1, 0)
            ret_material = material
        elif not isinstance(material, EchoSourceMaterial):
            raise TypeError(f"Attempted to add a material with type "
                                f"{type(material)}; must be an "
                                "EchoSourceMaterial or string when using "
                                f"dimensionality of {dimensionality}.")
        if dimensionality == "percent" or dimensionality == "fraction":
            self.materials.append((material, final_conc, "percent"))
        elif dimensionality == "volume":
            self.materials.append((material, final_conc, "volume"))
        else:
            raise ValueError(f"Attempted to add a material using dimensionality "
                             f"'{dimensionality}'. Dimensionality must be one of"
                             " 'concentration', 'volume', or 'percent'.")
        self.finalized = False

        return ret_material

    def current_vol(self):
        '''
        Calculates the total volume of all of the materials currently in the
        mix. Use this to figure out how much volume the mix contains, given
        what has already been added.
        '''
        return sum([vol for name, vol in self.recipe(finalize = False)])

    def fill_with(self, material):
        '''
        Fill all unfilled volume in the mix with some material (usually
        water). If another material was assigned to fill this, it will
        be overwritten by this call.

        Note that only one material can be used as a fill material at a time. If
        there was another fill material, it will be overwritten by this call.
        '''
        self.fill_material = material
        self.finalized = False

    def finalize(self):
        '''
        Finish up things like addition of fill materials. Also checks reaction
        for consistency, throwing a ValueError if the mix is overfilled or
        otherwise in obvious error, and raising a Warning if the reaction is
        underfull.
        '''
        #Check for overfilled wells
        if abs(self.get_volume() - self.current_vol()) > 1e-4:
            error_string = self.__class__.__name__
            error_string += " has %d nL volume but contains %.2f nL of " \
                            % (self.get_volume(), self.current_vol())
            error_string += "ingredients:"
            for material, material_vol in self.get_material_volumes():
                error_string += "\n\t%d nL of %s" % (material_vol, material)
            raise ValueError(error_string)

        if self.fill_material:
            fill_volume         = self.get_volume() - self.current_vol()
            fill_mat_final_conc = self.fill_material.nM * fill_volume \
                                  / self.get_volume()
            self.add_material(self.fill_material, fill_mat_final_conc)

        current_vol = self.current_vol()
        if current_vol - self.get_volume() > 1e-4:
            error_string = self.__class__.__name__
            error_string += " has %d nL volume but contains %.2f nL of " \
                            % (self.get_volume(), current_vol)
            error_string += "ingredients:"
            for material, material_vol in self.get_material_volumes():
                error_string += "\n\t%d nL of %s" % (material_vol, material)
            raise ValueError(error_string)

        if current_vol < self.get_volume():
            warn_string = self.__class__.__name__
            warn_string += "has %d nL volume but only contains %.2f nL of " \
                            % (self.get_volume(), current_vol)
            warn_string += "ingredients. Are you sure you want to underfill " \
                            + "this reaction?"
            for material, material_vol in self.get_material_volumes():
                warn_string += "\n\t%d nL of %s" % (material_vol, material)
            warnings.warn(warn_string, Warning)

        self.finalized = True

    def recipe(self, finalize = True):
        '''
        Generator returning descriptors of what goes in the mix.

        Arguments:
            finalize -- Iff True (default), and if the reaction hasn't already
                            been finalized, makes sure that it gets finalized.

        Yields -- pairs of the form (material, vol), where 'material' is an
                    EchoSourceMaterial and 'vol' is the volume of that material
                    to add to the reaction, in nL.
        '''
        # Make sure everything's ready to go and materials have been requested.
        if finalize and not self.finalized:
            self.finalize()

        for (material, vol) in self.get_material_volumes():
            name = str(material)
            yield (name, vol)

    def get_volume(self):
        '''
        The required total volume of the mix. Use this to know what the target
        volume of the mix should be, regardless of what's actually in it right
        now.
        '''
        return self.vol

    def get_material_volumes(self):
        '''
        Generator returning objects representing what goes in the mix, along
        with their volumes.
        '''
        if self.get_volume() == None:
            raise ValueError("self.vol is None for object "+repr(self))

        for material, final_conc, unit in self.materials:
            vol = self._convert_to_vol(material, final_conc, unit)
            yield material, vol

    def _convert_to_vol(self, material, final_conc, unit):
        '''
        Helper function for figuring out how much volume to use of a currently-
        stored material.
        '''
        if unit == "concentration":
            return final_conc * self.get_volume() / material.nM
        elif unit == "percent":
            return final_conc * self.get_volume()
        elif unit == "volume":
            return final_conc


class WellReaction(AbstractMixture):
    '''
    A reaction in a well on an Echo destination plate. Has a well, has
    volumes that are rounded to Echo-compatible numbers, and has a concept of
    echo-pipetted vs. hand-pipetted materials.
    '''
    def __init__(self, rxn_vol, well):
        super(WellReaction, self).__init__(rxn_vol)
        self.well = well
        self.hand_pipetted = dict()

    def fill_with(self, material, pipette_by_hand = False):
        '''
        Fill all unfilled volume in the mix with some material (usually
        water). If another material was assigned to fill this, it will
        be overwritten by this call.

        Note that only one material can be used as a fill material at a time. If
        there was another fill material, it will be overwritten by this call.
        '''
        self.fill_material_hand_pipetted = pipette_by_hand
        super(WellReaction, self).fill_with(material)

    def add_material(self, material, final_conc, dimensionality = "concentration",
                     pipette_by_hand = False):
        '''
        Add a material at a known final concentration. Final concentrations are
        assumed to use the same units as the material (usually nM). Does NOT
        check final reaction volume -- that will not be checked until
        finalize is called (happens automatically when recipe is
        called).

        Rounding to Echo-compatible volumes occurs at this step, unless the
        material is added by hand.
        '''
        self.hand_pipetted[material] = pipette_by_hand

        dimensionality = dimensionality.lower()
        if dimensionality == "concentration":
            target_vol  = self.vol * final_conc / material.nM
        elif dimensionality == "percent":
            target_vol = self.vol * final_conc
        elif dimensionality == "volume":
            target_vol = final_conc
        else:
            raise ValueError(f"Attempted to add a material to a WellReaction "
                              "using dimensionality "
                             f"'{dimensionality}'. dimensionality must be one "
                             "of 'concentration', 'volume', or 'percent'.")

        if pipette_by_hand:
            actual_conc = target_vol * material.nM / self.vol
        else:
            actual_vol  = echo_round(target_vol)
            actual_conc = actual_vol * material.nM / self.vol

        return super(WellReaction, self).add_material(material, actual_conc,
                                                      "concentration")

    def finalize(self):
        '''
        Checks reaction for consistency, throwing a ValueError if the reaction
        is overfilled or otherwise in obvious error, and raising a Warning if
        the reaction is underfull. Then, if everything checks out, volume is
        requested from the reaction's EchoSourceMaterials.
        '''
        current_vol = self.current_vol()

        if self.fill_material:
            fill_volume         = self.vol - current_vol
            fill_mat_final_conc = self.fill_material.nM * fill_volume \
                                  / self.vol
            self.add_material(self.fill_material, fill_mat_final_conc,
                              pipette_by_hand=self.fill_material_hand_pipetted)

        current_vol = self.current_vol()
        if current_vol > self.vol:
            error_string = "Reaction "
            if self.well:
                error_string += "in well %s " % self.well
            error_string += "has %d nL volume but contains %.2f nL of " \
                            % (self.vol, current_vol)
            error_string += "ingredients:"
            for material, material_vol in self.get_material_volumes():
                error_string += "\n\t%f nL of %s" % (material_vol, material)
            raise ValueError(error_string)
        if current_vol < self.vol:
            warn_string = "Reaction "
            if self.well:
                warn_string += "in well %s " % self.well
            warn_string += "has %d nL volume but only contains %.2f nL of " \
                            % (self.vol, current_vol)
            warn_string += "ingredients. Are you sure you want to underfill " \
                            + "this reaction?"
            warnings.warn(warn_string, Warning)

        for material, vol in self.get_material_volumes():
            material.request_material(self.well, vol,
                                      self.hand_pipetted[material])

        self.finalized = True


class MixtureMaterial(AbstractMixture, EchoSourceMaterial):
    '''
    Represents a mixture of materials that will be added to a source plate.
    I.e., it is an EchoSourceMaterial that can contain multiple
    components.

    MixtureMaterials have a rxn_vol (reaction volume), which is the size, in nL,
    of the reaction that MixtureMaterial is going into on the destination plate.
    Default 10000 (10 uL).

    A MixtureMaterial can EITHER contain materials with concentrations in
    dimensionality of concentration of the FINAL REACTION ("final_concentration")
    OR other dimensionality ("concentration", "volume", "fraction"). If final
    concentration dimensionality is used, then the MixtureMaterial should be
    added to reactions without a concentration (it will figure out its own
    volume). Otherwise, it should be added to reactions with a concentration
    (otherwise it can't figure out how much of each ingredient to use).
    '''
    def __init__(self, name, concentration = 1, units = None, vol = None,
                 well = None, length = 0, recipe_excess = 1.0, rxn_vol = 10000):
        AbstractMixture.__init__(self, vol = vol, well = well,
                                 recipe_excess = recipe_excess)
        EchoSourceMaterial.__init__(self, name, concentration = concentration,
                                    length = length, units = units)
        self.rxn_vol = rxn_vol
        self.mode = None
        self.plate = None
        self.name = name

    def fill_with(self, material):
        '''
        Fill all unfilled volume in the mix with some material (usually
        water). If another material was assigned to fill this, it will
        be overwritten by this call.

        Can't do this if using materials with dimensionality of final concentration,
        since the MixtureMaterial will calculate its volume based on its
        constituent materials.
        '''
        if self.mode == "final_concentration" \
           and not self.fill_material is None:
            raise ValueError(("Tried to fill a MixtureMaterial using dimensionality of ",
                             f"final concentration with material {material}."))
        self.fill_material = material
        self.finalized = False

    def add_material(self, material, final_conc,
                     dimensionality = "final_concentration"):
        '''
        Wraps AbstractMixture's add_material function, along with an additional
        (default) dimensionality option "final_concentration", which sets the final
        concentration in the reaction the MixtureMaterial is going into (NOT the
        final concentration within the mix itself).
        '''
        dimensionality = dimensionality.lower()

        if isinstance(material, str):
            material = EchoSourceMaterial(material, 1, 0)
        elif not isinstance(material, EchoSourceMaterial):
            raise TypeError(f"Attempted to add a material with type "
                                f"{type(material)}; must be an "
                                "EchoSourceMaterial or string when using "
                                f"dimensionality of {dimensionality}.")

        if dimensionality == "final_concentration":
            if self.mode is None:
                if not self.fill_material is None:
                    raise ValueError((f"Tried to add material {material} to a "
                                      "MixtureMaterial with a fill material"
                                      f" {self.fill_material}."))
                self.mode = "final_concentration"
            elif self.mode == "mix_concentration":
                raise ValueError((f"Tried to add a material {material} with "
                                  "dimensionality of final concentration to a "
                                  "MixtureMaterial that uses dimensionality of "
                                  "mix concentration."))
            self.materials.append((material, final_conc, "final_concentration"))
        else:
            if self.mode is None:
                self.mode = "mix_concentration"
            elif self.mode == "final_concentration":
                raise ValueError((f"Tried to add a material {material} with "
                                  "dimensionality of mix concentration to a "
                                  "MixtureMaterial that uses dimensionality of "
                                  "final concentration."))
            super().add_material(material, final_conc, dimensionality)

        return material

    def one_rxn_recipe(self, finalize = True):
        '''
        Iterator returning descriptors of what goes in a SINGLE REACTION.

        Arguments:
            finalize -- Iff True (default), and if the reaction hasn't already
                            been finalized, makes sure that it gets finalized.

        Yields -- pairs of the form (material, vol), where 'material' is an
                    EchoSourceMaterial and 'vol' is the volume of that material
                    to add to the reaction, in nL.
        '''
        if finalize and not self.finalized:
            self.finalize()

        if self.mode == "mix_concentration":
            for (material, mat_vol) in self.get_material_volumes():
                name = str(material)
                one_rxn_vol = self.rxn_vol * mat_vol / self.get_volume()
                yield (name, one_rxn_vol)
        elif self.mode == "final_concentration":
            for material, final_conc, dimensionality in self.materials:
                name = str(material)
                one_rxn_vol = self.rxn_vol * final_conc / material.nM
                yield (name, one_rxn_vol)
        elif len(self.materials) == 0:
            raise StopIteration
        else:
            raise ValueError(f"MixtureMaterial '{self}' has invalid unit_type"+\
                             f"'{self.mode}' (must be " + \
                              "'final_concentration' or 'mix_concentration').")

    def current_vol_per_rxn(self):
        '''
        Returns the current PER REACTION volume based on the materials
        currently in the master mix.
        '''
        return sum([vol for name, vol in self.one_rxn_recipe(finalize = False)])

    def text_recipe(self):
        '''
        Returns a string describing the mixture and its ingredients.
        '''
        ret_str = "\n%s:" % (self.name)
        ret_str += "\n\tStock Concentration: %d %s" % (self.concentration,
                                                       self.units)
        ret_str += "\n\tTotal Volume: %.2f" % (self.get_volume())
        lines = []
        for material, volume in self.recipe():
            lines.append("\n\t\t%0.2f uL %s" % (volume/1000, str(material)))
        return ret_str + "".join(lines)

    def get_volume(self):
        '''
        The required total volume of the mix. Use this to know what the target
        volume of the mix should be, regardless of what's actually in it right
        now.
        '''
        if self.total_volume_requested == 0:
            raise Exception(("Trying to get the volume of a MixtureMaterial ",
                             "object before requesting any material. ",
                             "MixtureMaterial objects should only be used as ",
                             "EchoSourceMaterials; if you want an abstract ",
                             "mixture, use the AbstractMixture class."))
        return self.total_volume_requested

    def get_material_volumes(self):
        '''
        Generator returning (material, volume) pairs of each material in the
        mix, and the volume of that material.
        '''
        if self.mode == "mix_concentration":
            for material, final_conc, unit in self.materials:
                vol = self._convert_to_vol(material, final_conc, unit)
                yield material, vol
        elif self.mode == "final_concentration":
            for material, one_rxn_vol in self.one_rxn_recipe(finalize = False):
                total_vol = self.get_volume() * one_rxn_vol / self.current_vol_per_rxn()
                yield material, total_vol

    def _convert_to_vol(self, material, final_conc, unit):
        '''
        Helper function for figuring out how much volume to use of a currently-
        stored material.
        '''
        if unit == "final_concentration":
            return (final_conc * self.rxn_vol / material.nM) \
                    * (self.get_volume() / self.current_vol_per_rxn())
        else:
            return super()._convert_to_vol(self, material, final_conc, unit)


class TXTLMasterMix(MixtureMaterial):
    '''
    A MixtureMaterial with built-in TX-TL extract and buffer.
    '''
    def __init__(self, extract_fraction = 0.33, mm_excess = 1.1,
                 rxn_vol = 10000, extract_per_aliquot = 30000,
                 buffer_per_aliquot = 37000, txtl_fraction = 0.75):
        '''
        extract_fraction: This is the fraction of the final mix made up of
                            TX-TL extract. Default 0.33 (lowest protein
                            concentration).
        mm_excess: The ratio of master-mix-to-make to total-mix-needed, i.e.,
                        mm_excess=1.1 => Make 10 percent excess, to account for
                        pipetting loss.
        rxn_vol: Total volume of a single reaction using this master mix, in nL
        extract_per_aliquot: Volume of TX-TL extract in one aliquot, in nL.
                                Default 30000.
        buffer_per_aliquot: Volume of TX-TL buffer in one aliquot, in nL.
                                Default 37000.
        txtl_fraction: Fraction of the total reaction allocated to
                        (extract + buffer)
        '''
        MixtureMaterial.__init__(self, name = "txtl_mm", concentration = 1,
                 vol = None, rxn_vol = rxn_vol, well = None, length = 0,
                 recipe_excess = mm_excess, units = "x")

        self.nM = 1.0   # 1x
        self.concentration = 1.0 #1x
        self.well     = "Master Mix"
        self.rxn_vol = rxn_vol
        self.vol     = rxn_vol

        self.mm_excess = mm_excess
        self.extract_fraction = extract_fraction
        self.extract_per_aliquot = extract_per_aliquot
        self.buffer_per_aliquot = buffer_per_aliquot
        self.txtl_fraction = txtl_fraction
        self.buffer_fraction = self.txtl_fraction - self.extract_fraction
        txtl_extract = EchoSourceMaterial("Extract", 1)
        txtl_buffer  = EchoSourceMaterial("Buffer",  1)
        self.add_material(txtl_extract, self.extract_fraction,
                          dimensionality = "final_concentration")
        self.add_material(txtl_buffer, self.buffer_fraction,
                          dimensionality = "final_concentration")

    def finalize(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            super(TXTLMasterMix, self).finalize()

    def text_recipe(self):
        ret_str  = "\nMaster Mix (including %d%% excess):"\
                        %((self.recipe_excess-1) * 100)
        ret_str += "\n\tTubes of extract needed: %d" % \
                        math.ceil(self.n_extract_aliquots())
        ret_str += "\n\tTubes of buffer needed: %d" % \
                        math.ceil(self.n_buffer_aliquots())
        recipe = super().text_recipe()
        ret_str += "\n" + "\n".join(recipe.split("\n")[1:])

        return ret_str

    def n_extract_aliquots(self, ):
        '''
        Returns the number of extract aliquots required to make this master mix.

        self.total_volume_requested should be set before calling this function.
        Throws an AttributeError otherwise. total_volume_requested is set during
        a call to request_picklist, when picks are finalized. If
        total_volume_requested is not set, will return 0.
        '''
        for material, conc, unit in self.materials:
            if material.name == "Extract":
                extract_vol = self.total_volume_requested * conc \
                                / material.nM / self.txtl_fraction
                return extract_vol / self.extract_per_aliquot \
                        * self.recipe_excess
        return 0

    def n_buffer_aliquots(self, ):
        '''
        Returns the number of buffer aliquots required to make this master mix.
        self.total_volume_requested must be set before calling this function.
        Throws an AttributeError otherwise. total_volume_requested is set during
        a call to request_picklist, when picks are finalized. If
        total_volume_requested is not set, will return 0.
        '''
        for material, conc, unit in self.materials:
            if material.name == "Buffer":
                buffer_vol = self.total_volume_requested * conc \
                                / material.nM / self.txtl_fraction
                return buffer_vol / self.buffer_per_aliquot \
                        * self.recipe_excess

        return 0


class SourcePlate():
    '''
    One Echo source plate. Responsible for allocating wells for
    EchoSourceMaterials.
    '''
    def __init__(self, filename = None, SPname = None, SPtype = None,
                 reuse_wells = False):
        '''
        SPname -- A(n arbitrary) string representing this plate. Default is
                    "Plate[1]"
        SPtype -- A string denoting the type of plate. Plate types are
                    restricted and determine the size of the plate. Default is
                    "384PP_AQ_BP"
        filename -- Name of a file holding a list of used wells for this plate.
                    If that file doesn't exist yet, it will be created and
                    populated when write_to_file is called.
        reuse_wells -- if False (default): well ingredients, concentrations and
                        volumes are not stored in the .dat file.
                       if True: well ingredients, concentrations, and volumes
                        are stored in the .dat file. When loaded again, these
                        wells are automatically re-used if they match user added
                        source materials.
        '''
        self.reuse_wells = reuse_wells
        self.materials_to_add = {} # stores materials added to wells for
                                   # instruction printing purposes
        if SPname == None:
            self.name = "Source[1]"
        else:
            name_is_number = True
            for c in SPname:
                if c not in string.digits:
                    name_is_number = False
                    break
            if name_is_number:
                self.name = "Source[%s]" % SPname
            else:
                self.name = SPname
        if SPtype == None:
            self.type = "384PP_AQ_BP"
            self.rows = 16
            self.cols = 24
        elif SPtype.startswith("384PP"):
            self.rows = 16
            self.cols = 24
            self.type = SPtype
        elif SPtype.startswith("6"):
            self.rows = 2
            self.cols = 3
            self.type = SPtype
        else:
            raise ValueError("'%s' is not a recognized plate type." % SPtype)
        self.wells_used   = {} #well(str)-->name, concentration, volume, date
        self.materials = {}#(well, name, concentration)-->[...(volume, date)...]
        self.wells_to_fill = {} #well --> (material, volume)
        self.current_row  = 0
        self.current_col  = 0

        self.used_well_file = filename
        if self.used_well_file:
            if os.path.isfile(filename):
                self.load_from_file(filename, reuse_wells = reuse_wells)

    def get_location(self, name, conc=None, i=0):
        """
        Get well location of material.

        Arguments:
            self: object
            name: name of material desired
            conc: select for concentration when available. Default is to select
                    the first entry.
        Returns:
            Well location, as a string
        Raises:
            AttributeError when named material is not found in source plate
        """
        def get_from_values(values):
            if type(values) is np.ndarray:
                n = min(len(values)-1, i)
                return values[n]
            return values

        names = self.plate.Name.values
        if name not in names:
            raise AttributeError(('No material named {} in source ' + \
                                  'plate').format(name))
        grouped_by_name = self.plate.groupby('Name')
        group = grouped_by_name.get_group(name)
        if conc is None:
            return get_from_values(group.Location.values)
        else:
            grouped_by_conc = group.groupby('Concentration')
            if conc in grouped_by_conc.groups:
                c_group = grouped_by_conc.get_group(conc)
                return get_from_values(c_group.Location.values)

    def load_from_file(self, filename, reuse_wells = False):
        '''
        Reads which wells have been used from a data file. The well-use file
        lists wells that have been used, with one line per well used.
        '''
        self.reuse_wells = reuse_wells
        with mt_open(filename, 'r') as infile:
            if not reuse_wells:
                well_ind = 0

            for line in infile:
                line = line.strip()
                if line == "":
                    continue
                elif line.lower()[:4]=="well" and not reuse_wells:
                    continue
                elif line.lower()[:4]=="well" and reuse_wells:
                    L = [l.strip() for l in line.lower().split(",")]
                    try:
                        well_ind = L.index("well")
                        name_ind = L.index("name")
                        conc_ind = L.index("concentration")
                        vol_ind = L.index("volume")
                        date_ind = L.index("date")
                    except ValueError:
                        raise ValueError("reuse_wells = True flag requires a "
                            ".dat file with a header line (csv format) "
                            "including the entries 'well', 'name', "
                            "'concentration', 'volume', 'date'")
                    continue

                L = line.split(",")
                well = L[well_ind]

                if self.reuse_wells:
                    name = L[name_ind]
                    try:
                        conc = float(L[conc_ind])
                    except ValueError:
                        conc = None
                    try:
                        vol = float(L[vol_ind])
                    except ValueError:
                        vol = None
                    try:
                        date = L[date_ind]
                    except ValueError:
                        date = None

                    self.wells_used[well] = (name, conc, vol, date)
                    if (name, conc) in self.materials:
                        self.materials[(name, conc)].append((well, vol, date))
                    else:
                        self.materials[(name, conc)] = [(well, vol, date)]
                else:
                    self.wells_used[well] = True

    def write_to_file(self):
        '''
        Writes which wells have been used to a data file. The well-use file
        lists wells that have been used, with one line per well used.

        Note that if no filename is given to the source plate, no well-use file
        will be written and this function will return silently.
        '''
        if not self.used_well_file:
            return

        with mt_open(self.used_well_file, 'w+') as outfile:
            if self.reuse_wells:
                outfile.write("well,name,concentration,volume,date\n")

            for well in self.wells_used:
                if self.reuse_wells:
                    (name, conc, vol, date) = self.wells_used[well]
                    outfile.write(f"{well},{name},{str(conc)},{vol},{date}\n")
                else:
                    outfile.write(well+"\n")

    def request_wells(self, n_wells):
        '''

        Called when an EchoSourceMaterial wants to get some wells. Returns a
        list of wells for the material, which are marked internally as used.

        Current logic is as follows: Scan from the top-left corner to the
        bottom-right corner, left to right, top to bottom, for a consecutive set
        of wells separated from used wells by a buffer well on the right and a
        buffer welll on the right. Assign the first block run across. If the
        number of wells requested is smaller than the number of wells per row,
        also require that the entire block be able to fit in one row.
        '''
        wells_used_array = np.zeros((self.rows, self.cols))
        for well in self.wells_used:
            row = well[0]
            row_ind = ord(row)-ord("A")
            col_ind = int(well[1:])-1
            wells_used_array[row_ind, col_ind] = 1

        if n_wells == 0:
            return []
        return_wells = np.empty(shape=(n_wells,), dtype=object)
        while True:
            # Check if this position will work
            if n_wells > self.cols or self.current_col + n_wells <= self.cols:
                flat_idx = self.current_row*self.cols + self.current_col
                if flat_idx + n_wells > self.rows * self.cols:
                    raise Exception("Source plate %s is out of available " + \
                                    "wells." % self.name)
                block = wells_used_array.ravel()[flat_idx:flat_idx + n_wells]
                # If it will, return that block and mark it used
                if True not in block:
                    block[:]    = True
                    return_list = np.empty(n_wells, dtype=object)
                    for i in range(flat_idx, flat_idx + n_wells):
                        row = int(i / self.cols)
                        col = i % self.cols
                        return_list[i-flat_idx] = string.ascii_uppercase[row] +\
                                         ("%02d" % (col+1))
                    # Mark wells as used.
                    self.increment_position(n_wells)
                    # Move on, leaving an extra well as a buffer (unless it
                    # just crossed to a new line)
                    if self.current_col != 0:
                        wells_used_array[self.current_row, self.current_col] = \
                                                                        True
                        self.increment_position(1)
                    return return_list
            # Otherwise, increment
            self.increment_position(1)

    def increment_position(self, n):
        '''
        DEPRICATED
        For internal use. Increments the plate's current position by n.
        '''
        if n < 0:
            raise ValueError("Can't increment a plate's position by negative "+\
                             "numbers. Good try, though.")
        self.current_col += n
        while self.current_col >= self.cols:
            self.current_col -= self.cols
            self.current_row += 1
        if self.current_row >= self.rows:
            raise Exception("Source plate %s is out of available wells." % \
                            self.name)

    def next_well(self, well):
        '''
        For internal use. Increments a well position by 1
        '''
        row = well[0]
        row_ind = ord(row)
        col = int(well[1:])

        if col <= self.cols:
            col+=1
            return row+str(col)
        elif row_ind < ord("A")+self.rows:
            row_ind+=1
            col = 1
            return chr(row_ind)+str(col)
        else:
            return None

    #Returns total available amount of a given material across all wells
    def get_available_material(self, material):
        name, conc = material.name, material.concentration
        tot_vol = 0
        if (name, conc) in self.materials:
            for (well, source_vol, date) in self.materials[name, conc]:
                tot_vol += source_vol-dead_volume
        return tot_vol

    def request_source_wells(self, material):
        """
        Returns the wells to fill with the given material as a list
        [(well, volume)]
        """
        #Material name and concentration
        name, conc = material.name, material.concentration
        #usable volume in a source well
        usable_volume  = max_volume - dead_volume
        #Total echo volume requested
        echo_volume = material.echo_volume_requested

        #Available material:
        tot_available_vol = self.get_available_material(material)
        #new material needed:
        volume_additional = echo_volume-tot_available_vol
        #Number of new source wells needed
        n_source_wells = math.ceil(float(volume_additional) /  usable_volume)

        #Fill Wells here
        wells_to_fill_list = self.request_wells(int(n_source_wells))
        wells_to_fill = []
        date = pydate.today().strftime("%d/%m/%Y")
        well_ind = 0

        while tot_available_vol < echo_volume:
            well = wells_to_fill_list[well_ind]
            if tot_available_vol+usable_volume >= echo_volume:
                fill_volume = echo_volume - tot_available_vol

                wells_to_fill += [(well, dead_volume+fill_volume)]

                tot_available_vol=echo_volume
                if (name, conc) in self.materials:
                    self.materials[name, conc] += \
                                    [(well, dead_volume+fill_volume, date)]
                else:
                    self.materials[name, conc] = \
                                        [(well, dead_volume+fill_volume, date)]
                self.wells_used[well] = (name, conc,
                                         dead_volume+fill_volume, date)
            elif tot_available_vol + usable_volume < echo_volume:

                wells_to_fill += [(well, max_volume)]

                tot_available_vol += usable_volume
                if (name, conc) in self.materials:
                    self.materials[name, conc] += [(well, max_volume, date)]
                else:
                    self.materials[name, conc] = [(well, max_volume, date)]
                self.wells_used[well] = (name, conc, max_volume, date)
            well_ind += 1
            material.total_volume_requested += dead_volume

        return wells_to_fill

    def request_picklist(self, material):
        '''
        Commit to wells, and return a list of finalized picks from this
        material.
        '''

        #Material name and concentration
        name, conc = material.name, material.concentration
        #usable volume in a source well
        usable_volume  = max_volume - dead_volume

        tot_vol_recieved = 0
        tot_vol_requested = 0
        picks = []
        #Iterate through picks (material --> Destinations)
        for pick in material.picklist:
            volume_requested = pick.volume
            tot_vol_requested += volume_requested
            #How much volume has been used so far for this pick
            volume_recieved = 0

            if not (name, conc) in self.materials:
                print(self.materials)
            for m_ind in range(len(self.materials[name, conc])):
                (well, source_vol, date) = self.materials[name, conc][m_ind]
                available_vol = source_vol-dead_volume
                #source well can fulfil the pick request
                if available_vol+volume_recieved >= volume_requested:
                    #Amount to send
                    volume_picked = volume_requested - volume_recieved
                    volume_recieved += volume_picked
                    #Update Source Plate Internal Data
                    self.materials[name, conc][m_ind] = \
                                          (well, source_vol-volume_picked, date)
                    self.wells_used[well] = (name, conc,
                                             source_vol-volume_picked, date)

                    #print(well,"-->",pick.destination_well, ":", volume_picked, "/", volume_recieved)
                    picks.append(Pick(material, well, pick.destination_well,
                                      volume_picked))
                    break
                #Additional source wells needed
                elif available_vol+volume_recieved < volume_requested \
                        and available_vol > 0:
                    #Amount to send
                    volume_picked = available_vol
                    volume_recieved += volume_picked

                    self.materials[name, conc][m_ind] = \
                                          (well, source_vol-volume_picked, date)
                    self.wells_used[well] = (name, conc,
                                             source_vol-volume_picked, date)

                    picks.append(Pick(material, well, pick.destination_well,
                                      volume_picked))
            tot_vol_recieved += volume_recieved

        return picks#, wells_to_fill

    def add_material_to_well(self, well, material, volume, update_date = False):
        """
        Adds a material to a well.
        Throws an error if the well is already full of a different material
        (or the same material at a different concentration).
        """


        name, conc = material.name, material.concentration
        if volume <= 0:
            warnings.warn(f"Adding volume={volume}<=0 of {name} to {well}. "
                          "Unecessary step omitted.")
            return

        date = pydate.today().strftime("%d/%m/%Y")

        if well not in self.wells_used:
            if volume > max_volume:
                raise ValueError(f"Cannot add volumes greater than "
                                 f"{max_volume} to source plate.")

            self.wells_used[well] = (name, conc, volume, date)
            if (name, conc) in self.materials:
                self.materials[(name, conc)]+=[(well, volume, date)]
            else:
                self.materials[(name, conc)] =[(well, volume, date)]

        elif self.wells_used[well][0] == name \
                    and self.wells_used[well][1] == conc:
            print("Refilling Well "+well+" with additional material "+name+".")
            old_vol = self.wells_used[well][2]
            old_date = self.wells_used[well][3]
            if volume+ old_vol> max_volume:
                raise ValueError("Cannot add volumes greater than "
                                 f"{max_volume-old_vol} to this well (which "
                                 f"already contains {old_vol}).")

            self.materials[(name, conc)].remove((well, old_vol, old_date))

            if update_date:
                self.wells_used[well] = (name, conc, volume+old_vol, date)
                self.materials[(name, conc)].append((well,
                                                     volume+old_vol, date))
            else:
                self.wells_used[well] = (name, conc, volume+old_vol, old_date)
                self.materials[(name, conc)].append((well, volume+old_vol,
                                                     old_date))

        else:
            raise ValueError(f"Well {well} already contains "
                             f"{self.wells_used[well][0]} at concentration "
                             f"{self.wells_used[well][1]}")

        material.plate = self

        if (name, conc) in self.materials_to_add:
            self.materials_to_add[(name, conc)]+=[(well, volume)]
        else:
            self.materials_to_add[(name, conc)] = [(well, volume)]



class DestinationPlate():
    '''
    One Echo destination plate. Responsible for allocating wells for
    reactions.
    '''
    def __init__(self, DPname = None, DPtype = None, filename = None):
        '''
        DPname -- A(n arbitrary) string representing this plate. Default is
                    "Destination[1]"
        DPtype -- A string denoting the type of plate. Plate types are
                    restricted and determine the size of the plate. Default is
                    "Nunc 384"
        filename -- Name of a file holding a list of used wells for this plate.
                    If that file doesn't exist yet, it will be created and
                    populated.
        '''
        if DPname == None:
            self.name = "Destination[1]"
        else:
            name_is_number = True
            for c in DPname:
                if c not in string.digits:
                    name_is_number = False
                    break
            if name_is_number:
                self.name = "Destination[%s]" % DPname
            else:
                self.name = DPname
        if DPtype == None or '384' in DPtype:
            self.type = "Nunc 384"
            self.rows = 16
            self.cols = 24
        elif '96' in DPtype:
            self.rows = 8
            self.cols = 12
            self.type = DPtype
        else:
            raise ValueError("'%s' is not a recognized plate type." % DPtype)

        # Construct use history object and helpers
        self.n_wells = self.rows * self.cols
        lets = string.ascii_uppercase[:self.rows]
        nums = range(1,1+self.cols)
        self.wells  = np.array(["{}{:02}".format(let, num) for let in lets \
                               for num in nums])
        self.indices = np.arange(self.n_wells)
        self.wells_used   = np.zeros((self.rows * self.cols,), dtype=bool)
        self.current_idx  = 0

        self.used_well_file = filename
        if self.used_well_file:
            if os.path.isfile(filename):
                self.load_from_file(filename)

    def make_picklist(self, source, rxns):
        n_rxns = len(rxns)
        d_wells = self.request_wells(n_rxns)
        picklist = []
        for rxn, well in zip(rxns, d_wells):
            for component in rxn:
                c_name = component[0]
                c_vol = component[1]
                if len(component)>2:
                    c_conc = component[2]
                else:
                    c_conc = None
                c_well = source.get_location(c_name, c_conc)
                picklist.append([c_well, c_name, well, c_vol, source.name,
                                 source.type, self.name])
        header = ["Source Well", "Sample Name", "Destination Well",
                  "Transfer Volume", "Source Plate Name", "Source Plate Type",
                  "Destination Plate Name"]
        return picklist, header

    def load_from_file(self, filename):
        '''
        Reads which wells have been used from a data file. The well-use file
        lists wells that have been used, with one line per well used.
        '''
        with mt_open(filename, 'r') as infile:
            for line in infile:
                line = line.strip()
                if line == "":
                    continue
                col_name = ""
                while line[-1] in string.digits:
                    col_name = line[-1] + col_name
                    line = line[:-1]
                col_num = int(col_name) - 1
                row_num = string.ascii_uppercase.find(line.upper())
                self.wells_used[row_num, col_num] = True

    def write_to_file(self):
        '''
        Writes which wells have been used to a data file. The well-use file
        lists wells that have been used, with one line per well used.
        '''

        if not self.used_well_file:
            return
        used_well_indices = self.wells_used.nonzero()
        used_well_rows = used_well_indices[0]
        used_well_cols = used_well_indices[1]
        with mt_open(self.used_well_file, 'w+') as outfile:
            for row_num, col_num in zip(used_well_rows, used_well_cols):
                row_string = string.ascii_uppercase[row_num]
                col_string = str(col_num + 1)
                outfile.write(row_string + col_string + "\n")

    def request_wells(self, n_wells):
        '''
        Called when an EchoSourceMaterial wants to get some wells. Returns a
        list of wells for the material, which are marked internally as used.

        Current logic is as follows: Scan from the top-left corner to the
        bottom-right corner, left to right, top to bottom, for a consecutive set
        of wells separated from used wells by a buffer well on the right and a
        buffer welll on the right. Assign the first block run across. If the
        number of wells requested is smaller than the number of wells per row,
        also require that the entire block be able to fit in one row.
        '''
        if n_wells == 0:
            return []
        if sum(self.wells_used) + n_wells > len(self.wells):
            raise Exception("Source plate %s is out of available wells." % \
                            self.name)
        unused_indices = self.indices[self.wells_used == False]
        return_indices = unused_indices[:n_wells]
        return_wells = self.wells[return_indices]
        self.wells_used[return_indices] = True
        return return_wells
