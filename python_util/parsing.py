from os.path import exists

class layer_section:
    def __init__(self, name):
        self.settings = dict()
        self.name = name

    def add_setting(self, key, value):
        self.settings[key] = value

    def get_configuration_value(self, key):
        parameter = None
        for string in self.settings.keys():
            if string == key:
                parameter = string
                break
        if parameter != None:
            return self.settings[parameter]

    def print_key_and_value(self):
        for key in self.settings.keys():
            print(key + " : " + self.settings[key])


def parsing_dataset(m_network):
    sections = parsing_config(m_network)
    for section in sections:
        if section.name == 'data':
            dataset = section.get_configuration_value("dataset") if section.get_configuration_value("dataset") != None else "imagenet"
        else :
            continue
    return dataset

def parsing_config(network):
    root_dir = '../'
    network_path = root_dir + 'configs/networks/' + network + '.cfg'
    
    # List for to store layers setting.
    sections = []
    file_exists = exists(network_path)
    if exists(network_path) == False :
        print('File does not exist at ' + str(network_path))
        exit()
    else :
        with open(network_path, 'r') as file:
            for line in file:
                string = line.strip().replace(" ", "")
                # Skipping black lines and comments
                if not string.strip() or string.startswith("#"):
                    continue
                # Discriminate layers from [section]
                if string.startswith("["):
                    section = layer_section(string[1:-1])
                    sections.append(section)
                # Layer configurations in the latest [section]
                else:
                    if string.find("=") == -1:
                        print("Error: Invalid config " + string)
                        exit()
                    key = string[:string.find("=")].lower()
                    value = string[string.find("=")+1:].lower()
                    sections[-1].add_setting(key, value)

    return sections
