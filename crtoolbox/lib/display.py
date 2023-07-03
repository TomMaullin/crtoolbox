import os
import numpy as np
from lib.fileio import read_image
from skimage import io
import plotly.graph_objects as go
import nibabel as nib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import colors

def display_volume(volume_file, mask = None, bg = None, mode='Sagittal', display='Fancy', slice=None):
    """
    """

    # If volume file is a string, read it in
    if isinstance(volume_file, str):
        # Read in volume
        volume = read_image(volume_file)
    else:
        # Set volume to volume_file
        volume = volume_file
        # Make sure it is a numeric numpy array
        volume = np.array(volume, dtype=np.float64)
    
    volume = volume.T
    
    # Remove any dimensions of size 1
    volume = np.squeeze(volume)

    # If we have a mask, load it in
    if mask is not None:

        # If mask is a string, read it in
        if isinstance(mask, str):
            # Read in mask
            vol = read_image(mask)
        else:
            # Set vol to mask
            vol = mask
            # Make sure it is a numeric numpy array
            vol = np.array(vol, dtype=np.int32)

        # Remove any dimensions of size 1
        vol = np.squeeze(vol)
        mask = vol.T

        # Mask out regions where we don't have data in the volume
        mask = mask * (np.abs(volume) > 1e-6)
    
    else:

        # Create mask
        mask = (np.abs(volume) > 1e-6)

    # If we have a background, load it in
    if bg is not None:

        # If bg is a string, read it in
        if isinstance(bg, str):              
            # Read in background
            vol = read_image(bg)
        else:
            # Set vol to bg
            vol = bg
            # Make sure it is a numeric numpy array
            vol = np.array(vol, dtype=np.float64)
        
        # Remove any dimensions of size 1
        vol = np.squeeze(vol)
        bg = vol.T

    else:

        # Create empty background
        bg = np.zeros(volume.shape)

    # Get the number of dimensions of the volume
    D = len(volume.shape)

    # Data type check
    mask = np.array(mask, dtype=np.int8)

    # If D is 2, reshape to a 3D volume and display in simple mode
    if D == 2:

        # Set mode to simple
        display = 'simple'
        
        # Set slice to 0
        slice = 0

        # Transpose volume
        volume = volume.T
        mask = mask.T
        bg = bg.T

        # Reshape to 3D
        volume = volume.reshape((1,) + volume.shape)
        mask = mask.reshape((1,) + mask.shape)
        bg = bg.reshape((1,) + bg.shape)

    else:

        # Check if we are in Sagital, Coronal or Axial view
        if mode.lower() == 'coronal':

            # Transpose the first two dimensions of the image
            volume = volume.transpose((1,0,2))
            mask = mask.transpose((1,0,2))
            bg = bg.transpose((1,0,2))

            # Flip the image along the second dimension
            volume = np.flip(volume,axis=1)
            mask = np.flip(mask,axis=1)
            bg = np.flip(bg,axis=1)

        elif mode.lower() == 'sagittal':

            # Transpose the first and last dimensions of the image
            volume = volume.transpose((2,0,1))
            mask = mask.transpose((2,0,1))
            bg = bg.transpose((2,0,1))

            # Flip the image along the second dimension
            volume = np.flip(volume,axis=1)
            mask = np.flip(mask,axis=1)
            bg = np.flip(bg,axis=1)

        elif mode.lower() == 'axial':

            # Flip the second dimension of the image
            volume = np.flip(volume,axis=1)
            mask = np.flip(mask,axis=1)
            bg = np.flip(bg,axis=1)

            # Flip the first dimension of the image
            volume = np.flip(volume,axis=0)
            mask = np.flip(mask,axis=0)
            bg = np.flip(bg,axis=0)

        else:

            # Raise error
            raise ValueError('Error: Invalid mode')

    # Check if the volume is entirely non-negative
    if np.all(volume >= 0):

        # Create a purely red color map
        cmap = 'Reds'

        # Colour for zero
        cmap_zero = cm.get_cmap('Reds')(0)

        # Multiply by 255 using a list comprehension
        cmap_zero = list((255*x for x in cmap_zero))

        # Remove transparency
        cmap_zero = cmap_zero

    # Check if the volume is entirely non-positive.
    elif np.all(volume <= 0):

        # Create a purely blue color map
        cmap = 'Blues'

        # Color for zero
        cmap_zero = cm.get_cmap('Blues')(1)

        # Multiply by 255 using a list comprehension
        cmap_zero = list((255*x for x in cmap_zero))

        # Remove transparency
        cmap_zero = cmap_zero

    # Otherwise, we have a mix of positive and negative values
    else:

        # Create a blue-red color map
        cmap = [[0, 'rgb(5,48,97)'],
                [0.1, 'rgb(33,102,172)'],
                [0.2, 'rgb(67,147,195)'],
                [0.3, 'rgb(146,197,222)'],
                [0.4, 'rgb(209,229,240)'],
                [0.5, 'rgb(247,247,247)'],
                [0.6, 'rgb(253,219,199)'],
                [0.7, 'rgb(244,165,130)'],
                [0.8, 'rgb(214,96,77)'],
                [0.9, 'rgb(178,24,43)'],
                [1.0, 'rgb(103,0,31)']]
                
        # Color for zero
        cmap_zero = (5,48,97)


    # Get dimensions
    nb_frames, c, h = volume.shape

    # Mask out zeros in volume
    volume[np.abs(mask) < 1e-6] = np.nan  

    # Add background
    volume[(np.abs(bg-1) < 1e-6)*np.isnan(volume)] = 0

    # Check if volume is all nan (i.e. no data)
    if np.all(np.isnan(volume)):
        # Default values
        vmin = 0
        vmax = 0
    else:
        # Compute global min and max
        vmin = np.nanmin(volume)  
        vmax = np.nanmax(volume)  

    # If the volume contains more than 2 unique values, set colorbar display to True
    if len(np.unique(volume)) > 3:
        colorbar = True
    else:
        colorbar = False

    # Check which display mode we are in
    if display.lower() == 'fancy':
            
        # Define frames
        frames = [
            go.Frame(
                data=[
                    go.Heatmap(
                        z=np.flipud(volume[nb_frames - 1 - k]),
                        zmin=vmin, zmax=vmax,
                        opacity=1,
                        showscale=colorbar,
                        colorscale=cmap,
                    ),
                    go.Heatmap(
                        z=[[vmin, vmax]],  # Invisible heatmap to force colorbar to appear on empty slices
                        opacity=0,
                        showscale=colorbar,
                        colorscale=cmap
                    )
                ],
                name=str(k))  # you need to name the frame for the animation to behave properly
            for k in range(nb_frames)
        ]

        # Add data to be displayed before animation starts
        fig = go.Figure(
            data=[
                go.Heatmap(
                    z=np.flipud(volume[nb_frames - 1]),
                    zmin=vmin, zmax=vmax,
                    opacity=1,
                    showscale=colorbar,
                    colorscale=cmap,
                ),
                go.Heatmap(
                    z=[[vmin, vmax]],  # Invisible heatmap  to force colorbar to appear on empty slices
                    opacity=0,
                    showscale=colorbar,
                    colorscale=cmap
                )
            ],
            frames=frames
        )

        # Add slider
        sliders = [
            {
                "pad": {"b": 10, "t": 60},
                "len": 0.9,
                "x": 0.1,
                "y": 0,
                "steps": [
                    {
                        "args": [[f.name], frame_args(0)],
                        "label": str(k),
                        "method": "animate",
                    }
                    for k, f in enumerate(fig.frames)
                ],
            }
        ]

        # If volume_file is a string, add title
        if isinstance(volume_file, str):
            # Make title using volume_file
            titlestr = 'Image: ' + volume_file.split(os.sep)[-1] + ' (' + mode + ' view)'
        else:
            # Make title
            titlestr = 'Image (' + mode + ' view)'

        # Layout
        fig.update_layout(
            title=titlestr,
            width=600,
            height=600,
            updatemenus=[
                {
                    "buttons": [
                        {
                            "args": [None, frame_args(50)],
                            "label": "&#9654;",  # play symbol
                            "method": "animate",
                        },
                        {
                            "args": [[None], frame_args(0)],
                            "label": "&#9724;",  # pause symbol
                            "method": "animate",
                        },
                    ],
                    "direction": "left",
                    "pad": {"r": 10, "t": 70},
                    "type": "buttons",
                    "x": 0.1,
                    "y": 0,
                }
            ],
            sliders=sliders,
            xaxis=dict(showticklabels=False, showgrid=False, zeroline=False),  # disable x-axis gridlines
            yaxis=dict(showticklabels=False, showgrid=False, zeroline=False),  # disable y-axis gridlines
            plot_bgcolor='rgb(0,0,0)',  # this sets the background color of the plot to black
        )

        # Display figure
        fig.show()

    # Check if we are in simple display mode
    elif display.lower() == 'simple':

        # Make sure we have a slice
        if slice is None:

            # Raise error
            raise ValueError('Error: Slice must be specified for simple display mode.')
        
        else:

            # Get slice
            img_slice = volume[nb_frames - 1 - slice,...]
            bg_slice = bg[nb_frames - 1 - slice,...]
            mask_slice = mask[nb_frames - 1 - slice,...]

            # create an empty image
            image = np.zeros((*img_slice.shape, 3))

            # Normalize the image
            if np.nanmax(img_slice) != np.nanmin(img_slice):
                img_slice = (img_slice - np.nanmin(img_slice)) / (np.nanmax(img_slice) - np.nanmin(img_slice))

            # Remove nans for now
            bg_slice[np.isnan(bg_slice)] = 0
            img_slice[np.isnan(img_slice)] = 0

            # assign the bg
            bg_slice = (np.abs(bg_slice - 1) < 1e-6)
            image[bg_slice, 0] = cmap_zero[0]/255 # red channel
            image[bg_slice, 1] = cmap_zero[1]/255 # green channel
            image[bg_slice, 2] = cmap_zero[2]/255 # blue channel

            # assign whats not in the bg
            image[~bg_slice, 0] = 0/255  # red channel
            image[~bg_slice, 1] = 0/255  # green channel
            image[~bg_slice, 2] = 0/255  # blue channel

            # Check if we have a color map already
            if isinstance(cmap, list):

                # Check if we need to divide by 255
                if not isinstance(cmap, str):
         
                    # Multiply by 255 using reformat_rgb function
                    cmap = reformat_rgb(cmap)

                # Convert the RGB values in list to a form matplotlib can use
                cmap = [colors.to_rgb(item[1]) for item in cmap]

                # Generate the colormap
                cmap = colors.LinearSegmentedColormap.from_list("my_colormap", cmap)

                # Apply colormap
                rgb_image = cmap(img_slice)

            else:

                # Apply colormap
                rgb_image = cm.get_cmap(cmap)(img_slice)

            # The returned array from colormap is a 3D array. The last dimension has 4 elements
            # (RGBA). We only want the first 3 (RGB)
            rgb_image = rgb_image[:, :, :3]

            # Assign the rgb into image
            mask_slice = (np.abs(mask_slice - 1) < 1e-6)
            image[mask_slice, 0] = rgb_image[mask_slice, 0]  # red channel
            image[mask_slice, 1] = rgb_image[mask_slice, 1]  # green channel
            image[mask_slice, 2] = rgb_image[mask_slice, 2]  # blue channel

            # display the image
            plt.imshow(image)
            plt.axis('off')  # hide the axis

            # Create a normalized color bar
            norm = colors.Normalize(vmin=vmin, vmax=vmax)

            # Create a ScalarMappable and initialize its data array
            sm = cm.ScalarMappable(norm=norm, cmap=cmap)
            sm.set_array(volume)

            # Add colorbar
            if colorbar:
                plt.colorbar(sm)
            
            # If its a 2D image, add a title
            if D == 2:

                # If volume_file is a string, add title
                if isinstance(volume_file, str):

                    # Make title using volume_file
                    titlestr = 'Image: ' + volume_file.split(os.sep)[-1]

                    # Add title
                    plt.title(titlestr)

            else:

                # If volume_file is a string, add title
                if isinstance(volume_file, str):
                    
                    # Make title using volume_file
                    titlestr = 'Image: ' + volume_file.split(os.sep)[-1] + ', Slice ' + str(slice) + ' (' + mode + ' view)'

                else:

                    # Make title
                    titlestr = 'Slice ' + str(slice) + ' (' + mode + ' view)'

                # Add title
                plt.title(titlestr)

            # show the image
            plt.show()

    else:

        # Raise error
        raise ValueError('Error: Invalid display mode. Options are "Fancy" or "Simple"')


def reformat_rgb(cmap):
    new_cmap = []
    for item in cmap:
        # Extract r, g, b values
        r, g, b = map(int, item[1][4:-1].split(',')) # this removes 'rgb(' and ')' and splits the values
        # Convert to [0, 1] range and append to new_cmap
        new_cmap.append([item[0], (r/255, g/255, b/255)])
    return new_cmap



def display_crs(fc, Upper_CR, Lower_CR, mask = None, mode='Sagittal', display='Fancy', slice=None):
    """
    Display confidence regions in 3D.

    Parameters:
    -----------
    fc : str
        Path to point estimate file.
    Upper_CR : str
        Path to upper confidence region file.
    Lower_CR : str
        Path to lower confidence region file.
    mask : str
        Path to mask file.
    mode : str
        View mode. Options are 'Sagittal', 'Coronal' or 'Axial'.

    Returns:
    --------
    None
    """

    # If upper confidence region is a string, read it in
    if isinstance(Upper_CR, str):
        # Read in volume
        vol = read_image(Upper_CR)
        volume_red = vol.T
    else:
        Upper_CR = 1*Upper_CR
        volume_red = Upper_CR.T

    # If lower confidence region is a string, read it in
    if isinstance(Lower_CR, str):
        # Read in volume
        vol = read_image(Lower_CR)
        volume_blue = vol.T
    else:
        Lower_CR = 1*Lower_CR
        volume_blue = Lower_CR.T

    # If point estimate is a string, read it in
    if isinstance(fc, str):
        # Read in volume
        vol = read_image(fc)
        volume_yellow = vol.T
    else:
        fc = 1*fc
        volume_yellow = fc.T

    # If we have a mask, load it in
    if mask is not None:

        # If mask is a string, read it in
        if isinstance(mask, str):
            # Read in mask
            vol = read_image(mask)
            volume_grey = vol.T
        else:
            mask = 1*mask
            volume_grey = mask.T
    
    else:

        # Create empty mask
        volume_grey = np.zeros(volume_red.shape)

    # Remove any dimensions of size 1
    volume_red = np.squeeze(volume_red)
    volume_blue = np.squeeze(volume_blue)
    volume_yellow = np.squeeze(volume_yellow)
    volume_grey = np.squeeze(volume_grey)

    # Get the number of dimensions of the volume
    D = len(volume_yellow.shape)

    # If D is 2, reshape to a 3D volume and display in simple mode
    if D == 2:

        # Set mode to simple
        display = 'simple'

        # Set slice to 0
        slice = 0
        
        # Transpose volume
        volume_grey = volume_grey.T
        volume_red = volume_red.T
        volume_yellow = volume_yellow.T
        volume_blue = volume_blue.T

        # Reshape to 3D
        volume_grey = volume_grey.reshape((1,) + volume_grey.shape)
        volume_red = volume_red.reshape((1,) + volume_red.shape)
        volume_yellow = volume_yellow.reshape((1,) + volume_yellow.shape)
        volume_blue = volume_blue.reshape((1,) + volume_blue.shape)

    else:

        # Check if we are in Sagital, Coronal or Axial view
        if mode.lower() == 'coronal':

            # Transpose the first two dimensions of the image
            volume_grey = volume_grey.transpose((1,0,2))
            volume_red = volume_red.transpose((1,0,2))
            volume_yellow = volume_yellow.transpose((1,0,2))
            volume_blue = volume_blue.transpose((1,0,2))

            # Flip the image along the second dimension
            volume_grey = np.flip(volume_grey,axis=1)
            volume_red = np.flip(volume_red,axis=1)
            volume_yellow = np.flip(volume_yellow,axis=1)   
            volume_blue = np.flip(volume_blue,axis=1)

        elif mode.lower() == 'sagittal':

            # Transpose the first and last dimensions of the image
            volume_grey = volume_grey.transpose((2,0,1))
            volume_red = volume_red.transpose((2,0,1))
            volume_yellow = volume_yellow.transpose((2,0,1))  
            volume_blue = volume_blue.transpose((2,0,1))

            # Flip the image along the second dimension
            volume_grey = np.flip(volume_grey,axis=1)
            volume_red = np.flip(volume_red,axis=1)
            volume_yellow = np.flip(volume_yellow,axis=1)
            volume_blue = np.flip(volume_blue,axis=1)

        elif mode.lower() == 'axial':

            # Flip the second dimension of the image
            volume_grey = np.flip(volume_grey,axis=1)
            volume_red = np.flip(volume_red,axis=1)
            volume_yellow = np.flip(volume_yellow,axis=1)
            volume_blue = np.flip(volume_blue,axis=1)

            # Flip the first dimension of the image
            volume_grey = np.flip(volume_grey,axis=0)
            volume_red = np.flip(volume_red,axis=0)
            volume_yellow = np.flip(volume_yellow,axis=0)
            volume_blue = np.flip(volume_blue,axis=0)

        else:

            # Raise error
            raise ValueError('Error: Invalid mode')

    # Get dimensions
    nb_frames, c, h = volume_grey.shape

    # Check which display mode we are in
    if display.lower() == 'fancy':

        # Define frames
        frames=[
            go.Frame(
                data=[
                    go.Heatmap(
                        z=np.flipud(volume_grey[nb_frames - 1 - k]),
                        zmin=0, zmax=1,
                        opacity=1,
                        showscale=False,
                        colorscale=[[0, 'rgba(0,0,0,1)'], [1, 'rgba(211, 211, 211, 1)']]),
                    go.Heatmap(
                        z=np.flipud(volume_blue[nb_frames - 1 - k]),
                        zmin=0, zmax=1,
                        opacity=1,
                        showscale=False,
                        colorscale=[[0, 'rgba(0,0,0,0)'], [1, 'blue']]),
                    go.Heatmap(
                        z=np.flipud(volume_yellow[nb_frames - 1 - k]),
                        zmin=0, zmax=1,
                        opacity=1,
                        showscale=False,
                        colorscale=[[0, 'rgba(0,0,0,0)'], [1, 'yellow']]),
                    go.Heatmap(
                        z=np.flipud(volume_red[nb_frames - 1 - k]),
                        zmin=0, zmax=1,
                        opacity=1,
                        showscale=False,
                        colorscale=[[0, 'rgba(0,0,0,0)'], [1, 'rgba(255, 0, 0, 1)']]),
                ],
                name=str(k))  # you need to name the frame for the animation to behave properly
            for k in range(nb_frames)]

        # Add data to be displayed before animation starts
        fig = go.Figure(
            data=[
                go.Heatmap(
                    z=np.flipud(volume_grey[nb_frames - 1]),
                    colorscale=[[0, 'rgba(0, 0, 0, 1)'], [1, 'rgba(211, 211, 211, 1)']],
                    zmin=0, zmax=1,
                    opacity=1,
                    showscale=False,
                ),
                go.Heatmap(
                    z=np.flipud(volume_blue[nb_frames - 1]),
                    colorscale=[[0, 'rgba(0,0,0,0)'], [1, 'blue']],
                    zmin=0, zmax=1,
                    opacity=1,
                    showscale=False,
                ),
                go.Heatmap(
                    z=np.flipud(volume_yellow[nb_frames - 1]),
                    colorscale=[[0, 'rgba(0,0,0,0)'], [1, 'yellow']],
                    zmin=0, zmax=1,
                    opacity=1,
                    showscale=False,
                ),
                go.Heatmap(
                    z=np.flipud(volume_red[nb_frames - 1]),
                    colorscale=[[0, 'rgba(0,0,0,0)'], [1, 'rgba(255, 0, 0, 1)']],
                    zmin=0, zmax=1,
                    opacity=1,
                    showscale=False,
                )
            ],
            frames=frames)

        # Add slider
        sliders = [
                    {
                        "pad": {"b": 10, "t": 60},
                        "len": 0.9,
                        "x": 0.1,
                        "y": 0,
                        "steps": [
                            {
                                "args": [[f.name], frame_args(0)],
                                "label": str(k),
                                "method": "animate",
                            }
                            for k, f in enumerate(fig.frames)
                        ],
                    }
                ]

        # Layout
        fig.update_layout(
                title='Confidence Regions (' + mode + ' view)',
                width=600,
                height=600,
                updatemenus = [
                    {
                        "buttons": [
                            {
                                "args": [None, frame_args(50)],
                                "label": "&#9654;", # play symbol
                                "method": "animate",
                            },
                            {
                                "args": [[None], frame_args(0)],
                                "label": "&#9724;", # pause symbol
                                "method": "animate",
                            },
                        ],
                        "direction": "left",
                        "pad": {"r": 10, "t": 70},
                        "type": "buttons",
                        "x": 0.1,
                        "y": 0,
                    }
                ],
                sliders=sliders,
                xaxis=dict(showticklabels=False),
                yaxis=dict(showticklabels=False),
        )

        # Display figure
        fig.show()
    
    # Check if we are in simple display mode
    elif display.lower() == 'simple':

        # Make sure we have a slice
        if slice is None:

            # Raise error
            raise ValueError('Error: Slice must be specified for simple display mode.')

        else:

            # Get slice
            slice_grey = volume_grey[nb_frames - 1 - slice,:,:]
            slice_red = volume_red[nb_frames - 1 - slice,:,:]
            slice_yellow = volume_yellow[nb_frames - 1 - slice,:,:]
            slice_blue = volume_blue[nb_frames - 1 - slice,:,:]

            # create an empty image
            image = np.zeros((*slice_grey.shape, 3))

            # assign the grey
            image[..., 0] = 211/255*slice_grey  # red channel
            image[..., 1] = 211/255*slice_grey  # green channel
            image[..., 2] = 211/255*slice_grey  # blue channel

            # assign the blue
            mask_blue = (np.abs(slice_blue - 1) < 1e-6)
            image[mask_blue, 0] = 0/255  # red channel
            image[mask_blue, 1] = 0/255  # green channel
            image[mask_blue, 2] = 255/255  # blue channel

            # assign the yellow
            mask_yellow = (np.abs(slice_yellow - 1) < 1e-6)
            image[mask_yellow, 0] = 255/255  # red channel
            image[mask_yellow, 1] = 255/255  # green channel
            image[mask_yellow, 2] = 0/255  # blue channel

            # assign the red
            mask_red = (np.abs(slice_red - 1) < 1e-6)
            image[mask_red, 0] = 255/255  # red channel
            image[mask_red, 1] = 0/255  # green channel
            image[mask_red, 2] = 0/255  # blue channel



            # display the image
            plt.imshow(image)
            plt.axis('off')  # hide the axis

            # If its a 2D image, add a title
            if D == 2:

                # Add title
                plt.title('Confidence Regions')

            else:

                # add a title
                plt.title('Confidence Regions, Slice ' + str(slice) + ' (' + mode + ' view)')

            # show the image
            plt.show()

    else:

        # Raise error
        raise ValueError('Error: Invalid display mode. Options are "Fancy" or "Simple"')
    

def frame_args(duration):
    """
    Define the frame arguments for the animation.

    Parameters:
    ----------
    duration : int

    Returns:
    --------
    dict : dictionary of frame arguments
    """

    return {
            "frame": {"duration": duration},
            "mode": "immediate",
            "fromcurrent": True,
            "transition": {"duration": duration, "easing": "linear"},
        }
