import pytest
import os
import sys
import shutil
import csv
import murraylab_tools.biotek as mt_biotek
import pandas as pd

data_folder       = "example_biotek_data"
supplement_folder = "example_biotek_supplements"
tidy_data_folder  = "tidy_data"

all_experiments = [f[:f.rfind(".")] for f in os.listdir(data_folder)]
# if not f.startswith(".") and os.path.isfile(os.path.join(data_folder, f))

@pytest.mark.parametrize('experiment_name', all_experiments)
def test_tidying(experiment_name):
    '''
    Tests the tidy_biotek_data function. This will be run on every data file in
    the data folder ("example_biotek_data"). Procedure for this test for each 
    data file is:

    1) Identify the supplement for that data file. The supplement file must be 
        in the supplement folder ("example_biotek_supplements") and must have 
        the same base name as the data file, plus an appended "_supplemental".
    2) Tidy the data, saving the result in the "tidy_data" folder. 
    3) Read out the newly-tidied data.
    4) Check the resulting DataFrame for consistency:
        4.1) The DataFrame should not be empty!
        4.2) The DataFrame's well list should be exactly the same as the list of 
                wells from the supplement file.
    '''
    # Convert to tidy format
    data_file_name  = os.path.join(data_folder, experiment_name + ".csv")
    supplement_name = os.path.join(supplement_folder, 
                                   experiment_name + "_supplemental.csv")
    try:
        mt_biotek.tidy_biotek_data(data_file_name, supplement_name)
    finally:
        tidy_file_name = experiment_name + "_tidy.csv"
        shutil.move(os.path.join(data_folder, tidy_file_name),
                    os.path.join(tidy_data_folder, tidy_file_name))

    # Check that the conversion worked.
    df = pd.read_csv(os.path.join(tidy_data_folder, tidy_file_name))
    assert len(df) > 0, "No data read for experiment" + experiment_name + "."

    well_list = []
    with open(supplement_name, 'r') as supp_file:
        reader = csv.reader(supp_file)
        next(reader)
        for line in reader:
            well_list.append(line[0])
    df_wells = df.Well.unique()
    assert [w in df_wells for w in well_list], \
            f"In experiment {experiment_name}, at least one well from the  " + \
            f" supplement doesn't appear in the corresponding tidied dataframe."
