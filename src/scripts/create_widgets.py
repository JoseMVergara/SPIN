"""
Author: @JoseMVergara
Description:
"""


import ipywidgets as widgets
from ipywidgets import interact, interactive, fixed, interact_manual, Button, Layout, jslink, IntText, GridBox, Label, Box, VBox
from IPython.display import display

def create_info_widget(tooltip):
    """
    Function to create helper widget
    info: help text
    """
    return Button(description='',
                        disabled=True,
                        button_style='', # 'success', 'info', 'warning', 'danger' or ''
                        tooltip=tooltip,
                        icon='info', # (FontAwesome names without the `fa-` prefix)
                        layout=Layout(height='25px', width='25px'))

def create_expanded_button(description, button_style):
    """
    Function to created expanded button widget
    input: button description, button style

    """
    return Button(description=description, button_style=button_style, layout=Layout(height='25px', width='auto'))

def create_int_text(value, description, step = 1):
    """
    Function to create int value text input widget
    inputs: int value, descripcion
    """
    return widgets.BoundedIntText(value=value,description=description,disabled=False, step=step,max=2000, min=0,layout=Layout(height='20px', width='50px'))

def create_float_text(value, description):
    """
    Function to create float value text input widget
    input: float value, description

    """
    return widgets.FloatText(value=value,description=description,disabled=False,layout=Layout(height='25px', width='70px'))

def create_float_slider(value, description):
    """
    Function to create a float slider widget
    input: float value, description
    """
    return widgets.FloatSlider(value=value,min=0,max=5,step=0.01,description=description,disabled=False,
                               continuous_update=False,orientation='horizontal',readout=True,readout_format='.1f',
                              layout=Layout(height='25px', width='70px'))

def create_text_input(value, description, width='150px'):
    """
    Function to create text input widget
    Input: text value, description
    """
    return widgets.Text(value=value, placeholder='Ingresa ',description=description,disabled=False,layout=Layout(height='25px', width=width))

def create_dropdown_input(options, value, description,width = '70px', disable=False):
    """
    Function to create dropdown input widget
    Input: options, default value, description  
    """
    return widgets.Dropdown(options=options,value=value,description=description,disabled=disable,layout=Layout(height='25px', width=width))

def create_text_area_input(value, placeholder, description):

    """
    Function to create text area input widget
    Input: Value input, placeholder, description
    """
    return widgets.Textarea(value=value,placeholder=placeholder,description=description,disabled=False,layout= Layout(height= '150px',width="600px"))

def create_combobox_input(options, placeholder, description):
    """
    Function to create combobox input widget
    Input: option list, placeholder, description
    """
    return widgets.Combobox(placeholder=placeholder,options=options,
                            description=description,ensure_option=True,disabled=False,layout=Layout(height='25px', width='70px'))

def create_checkbox_input(description, disable = False):
    """
    Function to create checkbox input widget
    Input: description
    """
    return widgets.Checkbox(value=False,description=description,disabled=disable,indent=False,layout=Layout(height='25px', width='70px'))

def create_toggle_input(options, description, style, tooltips):
    """
    Function to create toggle input widget
    Input: options list, description, style, tooltips list
    """
    return widgets.ToggleButtons(options=options,description=description,disabled=False,
                                 button_style=style,
                                 tooltips=tooltips,layout=Layout(height='25px', width='70px'))

def create_boolean_input(description, tooltip):
    """
    Function to create Boolean input widget
    Input: Description, tooltip
    """
    return widgets.ToggleButton(
    value=False,
    description=description,
    disabled=False,
    button_style='',
    tooltip=tooltip,layout=Layout(height='25px', width='150px'))

def HBox(list_):
    """
    Create Hbox Widget
    Input: List of  widgets
    """
    box_layout = widgets.Layout(display='flex',
                    flex_flow='row',
                    align_items='stretch',
                    border='doted',
                    width='90%')
    
    return widgets.HBox(list_, Layout = box_layout)