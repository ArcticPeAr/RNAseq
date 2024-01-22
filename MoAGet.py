 
import pubchempy as pcp

def get_cid_by_name(name):
    """Retrieve CID for a given molecule name using PubChemPy."""
    try:
        compounds = pcp.get_compounds(name, 'name', record_type='compound')
        # Usually, the first compound in the list is the most relevant, but be aware that multiple compounds can be returned
        return compounds[0].cid if compounds else None
    except Exception as e:
        print(f"Error retrieving CID for {name}: {e}")
        return None
    
molecule_name = "aspirin"
print(get_cid_by_name(molecule_name))