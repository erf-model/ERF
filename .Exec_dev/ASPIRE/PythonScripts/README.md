# ARM Data Retrieval & Plotting Tool

This repository provides instructions and scripts for downloading and visualizing atmospheric data from the Atmospheric Radiation Measurement (ARM) User Facility.

---

## 🛠️ Prerequisites & Setup

### 1. Create an ARM Account & Obtain a Token
To download data from the ARM Facility, you need a free user account and an API access token.

1. **Register for an account:** [ARM User Registration](https://adc.arm.gov/armuserreg/new)
2. **Log in to obtain your token:** [ARM SSO Login](https://sso.arm.gov/arm/login)

---

### 2. Install Dependencies

Install the `armlive_getfiles` utility directly from GitHub:

```bash
pip install git+[https://code.ornl.gov/ofg/armlive_getfiles.git](https://code.ornl.gov/ofg/armlive_getfiles.git)
```

```bash
getARMFiles -u <username>:<token> -ds <datastream> -s <start_date> -e <end_date>
```

```bash
python3 PlotARMData.py --arm-data-dir=<path-to-arm-data-dir>
```
