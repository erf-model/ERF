import cdsapi

# --- Read user values from file ---
def read_user_input(filename):
    inputs = {}
    with open(filename, 'r') as f:
        for line in f:
            if ':' not in line:
                continue
            key, value = line.strip().split(':', 1)
            key = key.strip()
            value = value.strip()
            if key == "area":
                inputs[key] = [float(x) for x in value.split(',')]
            else:
                inputs[key] = [value]
    return inputs

def download_era5_data(inputs):

    # Load values from user input file
    user_inputs = read_user_input(inputs)

    # Define dataset and static request fields
    dataset = "reanalysis-era5-pressure-levels"
    request = {
        "product_type": ["reanalysis"],
        "variable": [
            "geopotential",
            "relative_humidity",
            "specific_cloud_ice_water_content",
            "specific_cloud_liquid_water_content",
            "specific_humidity",
            "specific_rain_water_content",
            "specific_snow_water_content",
            "temperature",
            "u_component_of_wind",
            "v_component_of_wind",
            "vertical_velocity",
            "vorticity"
        ],
        "pressure_level": [
            "1", "2", "3", "5", "7", "10",
            "20", "30", "50", "70", "100", "125",
            "150", "175", "200", "225", "250", "300",
            "350", "400", "450", "500", "550", "600",
            "650", "700", "750", "775", "800", "825",
            "850", "875", "900", "925", "950", "975", "1000"
        ],
        "data_format": "grib",
        "download_format": "unarchived"
    }

    # Merge user input into request
    request.update(user_inputs)

    print("User inputs:", user_inputs)
    area = user_inputs.get("area", None)
    assert area and len(area) == 4, "'area' must be a list of four values"
    assert area[0] > area[2], "Latitude order invalid: North (1st) must be greater than South (3rd)"
    assert area[3] > area[1], "Longitude order invalid: East (4th) must be greater than West (2nd)"



    # Download data
    client = cdsapi.Client()
    filename = client.retrieve(dataset, request).download()
    print(f"Downloaded file: {filename}")
    return filename



