

# Instalar Gradio si no está instalado
!pip install gradio pandas chembl_webresource_client

# Importar bibliotecas necesarias
import pandas as pd
from chembl_webresource_client.new_client import new_client
import gradio as gr

# Función para realizar la búsqueda y filtrado en ChEMBL
def search_and_filter(target_name, sort_by_potency):
    # Crear un cliente para acceder a la API de ChEMBL
    chembl = new_client

    # Realizar la búsqueda del objetivo especificado
    target = chembl.target
    target_query = target.search(target_name)
    targets = pd.DataFrame.from_dict(target_query)

    if targets.empty:
        return "No se encontraron resultados para el objetivo especificado."

    # Obtener los datos de actividad biológica para el objetivo seleccionado
    activity = chembl.activity
    activity_query = activity.filter(target_chembl_id=targets['target_chembl_id'][0])
    activities = pd.DataFrame.from_dict(activity_query)

    if activities.empty:
        return "No se encontraron actividades para el objetivo especificado."

    # Ordenar los resultados por actividad EC50 si es necesario
    if sort_by_potency:
        activities.sort_values(by='standard_value', ascending=False, inplace=True)

    return activities


# Función de salida para mostrar los resultados
def output_results(target_name, sort_by_potency):
    activities = search_and_filter(target_name, sort_by_potency)
    if isinstance(activities, pd.DataFrame):
        return activities.head(10)
    else:
        return activities


# Lista de ejemplos predefinidos como listas de listas
examples_list = [["PAFR", False], ["AChE", True], ["TUBB", False]]

# Crear la interfaz de Gradio con los ejemplos corregidos
gr.Interface(fn=output_results,
             inputs=[target_name_input, sort_by_potency_input],
             outputs=gr.Dataframe(),
             examples=examples_list).launch()