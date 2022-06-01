import sys
import os
import json

def openJson(file):
    try:
        with open(file) as json_data:
            f = json.load(json_data)
    except IOError:
        print("Could not open the JSON config file: "+file)
        raise
    except ValueError:
        print("File format Error in the JSON config file: "+file)
        raise
    except Exception:
        print("Error while opening the JSON file: "+file)
        raise
    else:
        return f

json_file = sys.argv[1]
directory = json_file.rsplit("/",1)[0] #Remove .json from name
try:
    f = openJson(json_file)
except:
    exit(1)
if "bad" in json_file:
    print("Ignoring "+json_file)
    exit(0)
if "enabled" in f:
    if not f["enabled"]:
        print("Ignoring "+json_file)
        exit(0)
fields = ["authors", "name", "comparisons", "description", "executions", "postprocess"]
for field in fields:
    try:
        f[field]
    except:
        print("Missing field: "+field)
        exit(2)
if not f["name"] in directory:
    print("Incorrect name test: "+f["name"])
    exit(3)
if not os.path.isdir(directory):
    print("Directory does not exist: "+directory)
    exit(4)
dat = f["name"] + ".dat"
#if not os.path.isfile(dat):
#    print("Dat file does not exist: "+dat)
#    exit(5)
print(json_file+" format OK!")

