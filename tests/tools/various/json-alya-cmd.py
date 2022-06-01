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

def saveJson(d, name):
    """
        Save a json file
    """
    if(type(d) != type({})):
        logger.error("Date type is not a dict")
        raise TypeError("First parameter data type expected: dict " +
                        "but find: "+type(d))
    try:
        with open(name, 'w') as fp:
            json.dump(d, fp, indent=4, sort_keys=True)
    except IOError:
        logger.error("Could not open the JSON config file: "+name)
        raise
    except Exception:
        logger.error("Error while opening the JSON file: "+name)
        raise
    else:
        return 1

json_file = sys.argv[1]
f = openJson(json_file)
# Executions
name = f["name"]
a = []
status = False
for e in f["executions"]:
    if "commands" in e:
        if "alya" in e["commands"]:
            status = True
            cmd = e["commands"]["alya"]
            cmd = cmd.replace("mpirun", "[RUN]")
            cmd = cmd.replace("-np", "[NP]")
            cmd = cmd.replace("-n ", "[NP] ")
            cmd = cmd.replace(":", "[SEP]")
            cmd = cmd.replace(name, "[NAME]")
            e["commands"]["alya"] = cmd
        if "postprocess" in e["commands"]:
            status = True
            cmd = e["commands"]["postprocess"]
            cmd = cmd.replace(name, "[NAME]")
            e["commands"]["postprocess"] = cmd
    a.append(e)
f["executions"] = a
saveJson(f, json_file)
if status:
    print("OK")
