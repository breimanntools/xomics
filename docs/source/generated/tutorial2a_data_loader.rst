Data Loading Tutorial
=====================

This is a tutorial on loading of example proteomics datasets.

Various omics datasets can be loaded with ``xo.load_dataset``:

.. code:: ipython2

    import xomics as xo
    df_datasets = xo.load_dataset()
    df_datasets




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>Dataset</th>
          <th>Data Type</th>
          <th>Description</th>
          <th>Condition</th>
          <th>Quantification</th>
          <th>Reference</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>PROT_DEMYELINATION</td>
          <td>Proteomic</td>
          <td>Demylination recovery experiment in mouse</td>
          <td>4 timepoints in days [‘d00’, ‘d04’, ‘d10’, ‘d14’]</td>
          <td>LFQ</td>
          <td>Penktert21</td>
        </tr>
        <tr>
          <th>1</th>
          <td>LIPID_DEMYELINATION</td>
          <td>Lipidomics</td>
          <td>Demylination recovery experiment in mouse</td>
          <td>4 timepoints in days [‘d00’, ‘d04’, ‘d10’, ‘d14’]</td>
          <td>Standard</td>
          <td>Penktert21</td>
        </tr>
      </tbody>
    </table>
    </div>



.. code:: ipython2

    df_lfq = xo.load_dataset(name="LIPID_DEMYLINATION")
    df_lfq.head(5)





.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>timepoint</th>
          <th>sample</th>
          <th>ce_15:0;0</th>
          <th>ce_16:0;0</th>
          <th>ce_16:1;0</th>
          <th>ce_18:0;0</th>
          <th>ce_18:1;0</th>
          <th>ce_20:3;0</th>
          <th>ce_20:4;0</th>
          <th>ce_21:0;0</th>
          <th>...</th>
          <th>ps_16:1;0_22:5;0</th>
          <th>ps_18:2;0_20:4;0</th>
          <th>ps_16:1;0_22:6;0</th>
          <th>ps_18:0;0_21:0;0</th>
          <th>ps_17:0;0_22:1;0</th>
          <th>ps_17:1;0_22:0;0</th>
          <th>ps_18:0;0_21:1;0</th>
          <th>ps_18:1;0_21:0;0</th>
          <th>ps_19:0;0_20:1;0</th>
          <th>ps_17:1;0_22:1;0</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>0dpi</td>
          <td>0dpi-2</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>...</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>10.104308</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>3.241316</td>
          <td>NaN</td>
          <td>NaN</td>
        </tr>
        <tr>
          <th>1</th>
          <td>0dpi</td>
          <td>0dpi-3</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>...</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>0.054706</td>
          <td>NaN</td>
          <td>0.864512</td>
          <td>0.465531</td>
          <td>NaN</td>
        </tr>
        <tr>
          <th>2</th>
          <td>0dpi</td>
          <td>0dpi-4</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>...</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>5.341029</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>1.135322</td>
          <td>0.576992</td>
          <td>NaN</td>
        </tr>
        <tr>
          <th>3</th>
          <td>0dpi</td>
          <td>0dpi-5</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>...</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>8.250989</td>
          <td>0.12316</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>1.687771</td>
          <td>0.788158</td>
          <td>0.070025</td>
        </tr>
        <tr>
          <th>4</th>
          <td>3dpi</td>
          <td>3dpi-1</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>...</td>
          <td>0.809395</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>0.046467</td>
          <td>NaN</td>
          <td>0.689641</td>
          <td>0.286349</td>
          <td>NaN</td>
        </tr>
      </tbody>
    </table>
    <p>5 rows × 1024 columns</p>
    </div>



