# you need a config file as such : "example_cmapPy_config_file.cfg" but renamed to .cmapPy.cfg
import cmapPy
import pandas as pd

from cmapPy.clue_api_client.clue_api_client import ClueApiClient

ClueApiClient.run_count_query("down-reg_GO:1902430.xlsx",  )
