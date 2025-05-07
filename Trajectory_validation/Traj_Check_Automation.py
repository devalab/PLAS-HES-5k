#!/usr/bin/env python
from prody import *
import numpy as np
import os
import argparse
import sys
import subprocess
import concurrent.futures
from os.path import exists
from colorama import Fore, Style
from concurrent.futures import ThreadPoolExecutor
from google.oauth2.service_account import Credentials
from googleapiclient.discovery import build


# Variables
username = 'Usename'
ada_user='ada_username@ada'
absolute_path='path/to/files/for/checking'
SERVICE_ACCOUNT_FILE = 'path/to/api_json/file'


# Constants
SPREADSHEET_ID = 'SPREADSHEET_ID'
SHEET_NAME = 'HES'
SCOPES = ['https://www.googleapis.com/auth/spreadsheets']


BASH = "/bin/bash"
ZSH = "/bin/zsh"

# Authenticate Google Sheets API
credentials = Credentials.from_service_account_file(SERVICE_ACCOUNT_FILE, scopes=SCOPES)
service = build('sheets', 'v4', credentials=credentials)
sheet = service.spreadsheets()

# convert column index to letter
def col_index_to_letter(index):
    letter = ''
    while index > 0:
        index, remainder = divmod(index - 1, 26)
        letter = chr(65 + remainder) + letter
    return letter
    
# update the Google Sheet
def update_sheet(pdb_id, run_number, status):
    # D=4, E=5, F=6, G=7, H=8
    column_letter = col_index_to_letter(run_number + 5)
    sheet_range = f'{SHEET_NAME}!{column_letter}:{column_letter}'

    result = sheet.values().get(spreadsheetId=SPREADSHEET_ID, range=SHEET_NAME).execute()
    values = result.get('values', [])

    # Find row with matching PDB ID and check username in the third column
    row_index = next((i for i, row in enumerate(values) if row[1] == pdb_id and row[2] == username), None)
    if row_index is not None:
        update_range = f"{SHEET_NAME}!{column_letter}{row_index + 1}"
        body = {'values': [[status]]}
        request = sheet.values().update(spreadsheetId=SPREADSHEET_ID, range=update_range,
                                        valueInputOption='USER_ENTERED', body=body)
        response = request.execute()
        print(f"Updated PDB ID {pdb_id} at run {run_number} with status '{status}' for user '{ada_user}'.")
    else:
        print(f"PDB ID {pdb_id} not found in the sheet for user '{ada_user}' or the username does not match.")


#update IP column
def update_ip(pdb_id, run_number, status):
    # D=4, E=5, F=6, G=7, H=8
    column_letter = col_index_to_letter(run_number + 10)
    sheet_range = f'{SHEET_NAME}!{column_letter}:{column_letter}'

    result = sheet.values().get(spreadsheetId=SPREADSHEET_ID, range=SHEET_NAME).execute()
    values = result.get('values', [])

    # Find row with matching PDB ID and check username in the third column
    row_index = next((i for i, row in enumerate(values) if row[1] == pdb_id and row[2] == username), None)
    if row_index is not None:
        update_range = f"{SHEET_NAME}!{column_letter}{row_index + 1}"
        body = {'values': [[status]]}
        request = sheet.values().update(spreadsheetId=SPREADSHEET_ID, range=update_range,
                                        valueInputOption='USER_ENTERED', body=body)
        response = request.execute()
        print(f"Updated PDB ID {pdb_id} at run {run_number} with status '{status}' for user '{ada_user}'.")
    else:
        print(f"PDB ID {pdb_id} not found in the sheet for user '{ada_user}' or the username does not match.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run scripts")
    parser.add_argument("index", type=str, nargs="?")
    parser.add_argument("pdbid", type=str, nargs="?")
    parser.add_argument("-r", "--read_input", action="store_true")

    if parser.parse_args().read_input:
        input_list = f"{absolute_path}/input.txt"
    elif parser.parse_args().index and parser.parse_args().pdbid:
        args = parser.parse_args()
        index = args.index
        pdbid = args.pdbid
        input_list = f"{absolute_path}/{index}.{pdbid}.txt"
        with open(input_list, "w", encoding="utf-8") as f:
            f.write(f"{index} {pdbid}")
    else:
        print("Please enter the index and pdbid, or provide an input file")
        sys.exit()


    with open(input_list, "r", encoding="utf-8") as id_list:
        for line in id_list:
            if not line.strip():
                continue
            index, pdbid = line.strip().split()[0],line.strip().split()[1]

            if not index or not pdbid:
                continue
            input_path = f"{absolute_path}/{index}.{pdbid}"

            print('###############################################')
            message = 'Processing index {} and pdbid {}'.format(index,pdbid)
            highlighted_message = "{}{}{}".format(Fore.GREEN, message, Style.RESET_ALL)
            print(highlighted_message)

            
            os.chdir(f"{input_path}")
            
            message = 'Processing index {} and pdbid {} initial structure'.format(index,pdbid)
            highlighted_message = "{}{}{}".format(Fore.YELLOW, message, Style.RESET_ALL)
            print(highlighted_message)
           
                
            for i in range(5):
                
                message = 'Processing index {} and pdbid {} ---------> run {}'.format(index,pdbid,i)
                highlighted_message = "{}{}{}".format(Fore.GREEN, message, Style.RESET_ALL)
                print(highlighted_message)
                
                                
                os.chdir(f"{absolute_path}/{index}.{pdbid}/run{i}")
                                                #VMD Startup Command
                command = f"""/usr/local/bin/vmd -e ../../state_1.vmd"""
                subprocess.run(command, shell=True, check=True, executable=BASH)
                
                                # Prompt user for input
                user_input = input(f"Is the initial I/P run {i} for PDB ID {pdbid} correct? (y/n): ").strip().lower().format(Fore.GREEN, input, Style.RESET_ALL)
                status = 'Yes' if user_input == 'y' else 'No'
                # Update Google Sheet
                update_ip(pdbid, i, status)
                
                
                                                #VMD Startup Command
                command = f"""/usr/local/bin/vmd -e ../../state.vmd"""
                subprocess.run(command, shell=True, check=True, executable=BASH)
                command = f"""cp {absolute_path}/distance_1.py {absolute_path}/{index}.{pdbid}/run{i}/"""
                subprocess.run(command, shell=True, check=True, executable=BASH)
                command = f"""cp {input_path}/complex_solvated.prmtop {absolute_path}/{index}.{pdbid}/run{i}/"""
                subprocess.run(command, shell=True, check=True, executable=BASH)
                command = f"""python distance_1.py"""
                subprocess.run(command, shell=True, check=True, executable=BASH)
                
                # Prompt user for input
                user_input = input(f"Is the run {i} for PDB ID {pdbid} correct? (y/n): ").strip().lower().format(Fore.GREEN, input, Style.RESET_ALL)
                status = 'Yes' if user_input == 'y' else 'No'
                # Update Google Sheet
                update_sheet(pdbid, i, status)
                

            os.chdir(f"{absolute_path}")
            os.system('rm -rf {}.{}'.format(index,pdbid))


