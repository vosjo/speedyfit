import pytest

import pandas as pd

from speedyfit.speedyfit_batch import process_objects

class TestBatch:

    def test_process_objects(self):

        object_data = pd.DataFrame(data={
            'name':['GALEXJ153213.7+010451', 'HD42709', 'Gaia DR2 1209876314302933632'],
            'info':['', 'intersting star', 'bla']
        })

        object_list = process_objects(object_data)

        assert len(object_list) == 3
        assert object_list[2] == 'Gaia_DR2_1209876314302933632'

        object_data = pd.DataFrame(data={
            'ra': [233.057210996, 10.7429509085],
            'dec': [1.0807152108, -38.1270271017]
        })

        object_list = process_objects(object_data)

        assert len(object_list) == 2
        assert object_list[0] == 'J153213.73+010450.57'
        assert object_list[1] == 'J004258.31-380737.30'


        object_data = pd.DataFrame(data={
            'ra': ['15:32:13.73', '00:42:58.31'],
            'dec': ['+01:04:50.57', '-38:07:37.30']
        })

        object_list = process_objects(object_data)

        assert len(object_list) == 2
        assert object_list[0] == 'J153213.73+010450.57'
        assert object_list[1] == 'J004258.31-380737.30'

        object_data = pd.DataFrame(data={
            'ra': ['15 32 13.73', '00 42 58.31', '15 32 13.73'],
            'dec': ['+01 04 50.57', '-38 07 37.30', '01 04 50.57']
        })

        object_list = process_objects(object_data)

        assert len(object_list) == 3
        assert object_list[0] == 'J153213.73+010450.57'
        assert object_list[1] == 'J004258.31-380737.30'
        assert object_list[2] == 'J153213.73+010450.57'
