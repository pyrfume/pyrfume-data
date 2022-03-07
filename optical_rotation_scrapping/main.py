import re

def extract_oprvals(list_of_vals):
    """Takes in a dictionary of cids of molecules matched with strings that have the optical
    rotation value(s) for each molecule. Returns a list of optical rotation values for each
    CID in the dictionary. Some molecules have multiple strings and/or optical rotation values 
    given. In these case, the optical rotation value that was extracted at a temperature closest
    to 20 degrees (room temperature) was taken. If temperature values were not given then any 
    optical rotation value was choosen.

    Args:
        dict_of_vals ([dict]): Keys are cids and values are strings

    Returns:
        [list]: optical rotation values for each molecule
    """
    orv = []
    for val_string in list_of_vals:
        opr_values_found = {}
        for s in val_string:
            val = "".join(s.split())
            if '°C/' in val:
                out = list(parse_temperature_values(val, '°C/'))
                opr_values_found[out[0]] = out[1]
            elif '°C' in val:
                out = list(parse_temperature_values(val, '°C'))
                opr_values_found[out[0]] = out[1]
            else:
                opr_values_found["X"] = parse_opr_values(val)
        
        orv.append(get_opr_from_dict(opr_values_found))     
    return orv


def get_opr_from_dict(dict_of_opr_vals):
    """Takes in a dictionary where the keys are temperatures and values are optical rotation 
    values. The dictionary is for all the temperatures and optical rotation values extracted 
    for one molecule. This function determines which of the values in the dictionary to keep.

    Args:
        dict_of_opr_vals ([dict]): Keys are temperature and values are optical rotation vals.

    Returns:
        [String]: Final optical rotation value for a molecule
    """
    if len(dict_of_opr_vals) > 0:
        dict_keys = list(dict_of_opr_vals.keys())
        if dict_keys.count("") == len(dict_keys):
            return dict_of_opr_vals[""]
        if "" in dict_keys:
            dict_keys.remove("")
        if dict_keys.count("X") == len(dict_keys):
            return dict_of_opr_vals["X"]
        else:
            try:
                dict_keys.remove("X")
            except:
                pass
            return dict_of_opr_vals[dict_keys[abs_distance(dict_keys)]]
    else:
        return dict_of_opr_vals[0]
        

def parse_temperature_values(val_given, temp_symbol):
    """Given a string with temperature and optical rotation values for one molecule. Given a symbol
    that's used to parse through a find temperature values. Finds all the temperature values in one
    string. Returns the temperature value that closest to 20 degrees along with its respective 
    optical rotation value.

    Args:
        val_given ([str]): string of temperature and optical rotation values
        temp_symbol ([str]): Specific substring symbol to find temperature values

    Returns:
        [tuple]: temperature value, optical rotatation value associated w that temperature

    Notes:
        There should not be a case where single_temp_val or temp_val stay none. The reason why this 
        happens and there is code to account for it is because parsing messes up if the substring given is
        "°C" followed immediately after by "(" 
        So code needs to be fixed to get parsing to work for when "°C(" comes up
        Because of this error, None values show up in some in the end opr list
    """
    temperature_values = []
    deg_values = []

    if val_given.count(temp_symbol) > 1:
        while(1):
            if temp_symbol in val_given:
                temp_val = None
                try:
                    temp_val = re.search(fr"(\d+)\d{{0}}".format(temp_symbol), val_given).group()
                    temperature_values.append(temp_val[0:temp_val.index(temp_symbol)])
                    deg_values.append(parse_opr_values(val_given))
                except:
                    pass
                val_given = val_given[val_given.index(temp_symbol)+3: -1]
            else:
                break
    else:
        single_temp_val = None
        try:
            stv = re.search(fr"(\d+)\d{{0}}".format(temp_symbol), val_given).group()
            single_temp_val = stv[0:stv.index(temp_symbol)]
        except:
            pass

        if single_temp_val is None:
            return "", parse_opr_values(val_given)
        else: 
            return single_temp_val, parse_opr_values(val_given)

    index = abs_distance(temperature_values)
    return temperature_values[index], deg_values[index]


def parse_opr_values(val):
    """Takes in a string. Parses the string for an optical rotation value.
    Returns. optical rotation value otherwise nothing.    

    Args:
        val ([str]): string associated with a molecule that contains optical
        rotation values

    Returns:
        [str]: optical rotation value, or returns nothing
    """
    search_strings = ['\+\d+\.\d+', '\+\d+', '\-\d+\.\d+', '\-\d+']
    for s in search_strings:
        try:
            return_value = re.search(s, val)
            return return_value.group()
        except:
            pass


def abs_distance(temp_vals):
    """Takes a list of numbers and returns the idex of the number in the list
    that has the smallest distance to 20.

    Args:
        temp_vals ([temp_vals]): List of numbers

    Returns:
        [Int]: Index of the number in the list that is closest to 20

    Note:
        If multiple numbers in the list have the same distance to 20. The first 
        value that is reached during iteration is kept.
    """
    smallest_dist = float('inf')
    index_of_smallest = -1
    for i in range(0, len(temp_vals)):
        val = temp_vals[i]
        try:
            value = float(val)
        except:
            value = float(val[:-3])
        if abs(value - 20) < smallest_dist:
            smallest_dist = abs(value - 20)
            index_of_smallest = i
    return index_of_smallest


def cross_checks(human_nums, output_nums):
    """Recieves two array of numbers. Compares the two arrays and returns how many
    mismatches are between the two arrays. This function checks to see how identical
    the two arrays are.

    Args:
        human_nums ([list]): Human read optical rotation values
        output_nums ([list]): Code computed optical rotation values

    Returns:
        [Int]: How many values these two arrays differ by
    """
    count_diffs = 0
    for i in range(0, len(human_nums)):
        if human_nums[i] != output_nums[i]:
            count_diffs += 1
            print("This diff happens at ", i)
            print(human_nums[i])
            print(output_nums[i])
    return count_diffs