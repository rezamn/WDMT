from myParser import parse_my_line
from loadData import load_data, list_stream
from processData import clean_stream, station_list

args = parse_my_line()
print(args)
data = load_data(args)
data = clean_stream(args, data)
# list_stream(data)
# data.merge()
# gap_data = data.get_gaps()
list_stream(data)