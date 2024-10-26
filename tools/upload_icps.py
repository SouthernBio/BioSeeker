# upload_icps.py
import pandas as pd
import os

def upload_icps(file_path: str):
    """Function to upload ICPs from a CSV file and store it in the appropriate directory.

    Args:
        file_path (str): Path of the file to upload.

    Returns:
        str: Confirmation message.
    """
    if not os.path.exists(file_path) or not file_path.endswith('.csv'):
        return "The file does not exist or is not a valid CSV file."

    try:
        df = pd.read_csv(file_path)
        
        # Check if the file contains the necessary columns (adapt as needed)
        required_columns = ['Codon', 'ConservationRate']  # Example of expected columns
        if not all(col in df.columns for col in required_columns):
            return "The CSV file does not contain all required columns."

        # Save the file to a specific directory
        output_dir = 'user_uploaded_icps'
        os.makedirs(output_dir, exist_ok=True)
        df.to_csv(os.path.join(output_dir, 'uploaded_icps.csv'), index=False)

        return "ICPs file uploaded successfully."
    except Exception as e:
        return f"Error while uploading the ICPs file: {str(e)}"
