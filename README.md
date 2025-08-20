A simple Streamlit app showing how various parameters impact the transfer function of an accretion disk. The wavelength dependent transfer function represents the kernel between the driving variability and the response at specified wavelengths. Various sliders control the typical accretion disk parameters and our viewing conditions. Applying the "GR" toggle will apply relativistic corrections via pre-calculated raytraces computed using Sim5. 

[![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](<insert streamlit link here>)

### How to run it on your own machine

1. Install the requirements

   ```
   $ pip install -r requirements.txt
   ```

2. Run the app

   ```
   $ streamlit run streamlit_app.py
   ```

### How to expand this script to include more inclination angles (requires ~2.4 GiB storage on local device for pre-computed GR ray trace files)

1. Clone the github repository `$ git clone https://github.com/Henry-Best-01/visualize_disk_transfer_function.git`

2. Navigate to the cloned repository `$ cd visualize_disk_transfer_function`

3. Download `https://drive.google.com/file/d/1WGuzCXMKXwSawcfF8LE1gQTc042pzz_2/view?usp=share_link`

3. Unzip the GR ray trace data

4. In "streamlit_app.py", make the following changes to lines:
   
   15: path_to_raytraces = "<Path_to_your_downloaded_ray_traces>"
   
   24: inclination = left_col.slider("inclination angle [deg.]", min_value=0, max_value=89, step=1, value=20)

5. Follow instructions for "How to run it on your own machine" above.

