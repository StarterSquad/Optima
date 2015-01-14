#!/bin/env python
# -*- coding: utf-8 -*-
from optima_test_base import OptimaTestCase
import unittest
import json
from api import app
from flask import helpers

class ProjectTestCase(OptimaTestCase):
    """
    Test class for the project blueprint covering all /api/project endpoints.

    """

    def setUp(self):
        super(ProjectTestCase, self).setUp()
        response = self.create_user()
        response = self.login()

    def test_create_project(self):
        response = self.client.post('/api/project/create/test', data = '{}')
        self.assertEqual(response.status_code, 200)

    def test_retrieve_project_info_fails(self):
        headers = [('project', 'test')]
        response = self.client.get('/api/project/info', headers=headers)
        self.assertEqual(response.status_code, 500)
        self.assertEqual(json.loads(response.data), { "status": "NOK" })

    def test_retrieve_project_info(self):
        self.create_project('test')

        headers = [('project', 'test')]
        response = self.client.get('/api/project/info', headers=headers)
        self.assertEqual(response.status_code, 200)
        project_data = json.loads(response.data)
        self.assertEqual(project_data['name'], 'test')
        self.assertEqual(project_data['status'], 'OK')

    def test_retrieve_project_list(self):
        self.create_project('test2')

        response = self.client.get('/api/project/list')
        self.assertEqual(response.status_code, 200)
        projects_data = json.loads(response.data)
        self.assertEqual(projects_data['projects'][0]['name'], 'test2')
        self.assertEqual(projects_data['projects'][0]['status'], 'OK')

    def test_project_params(self):
        from sim.parameters import parameter_name
        response = self.client.get('/api/project/params')
        print(response)
        self.assertEqual(response.status_code, 200)
        params = json.loads(response.data)['params']
        self.assertTrue(len(params)>0)
        self.assertTrue(set(params[0].keys())== \
            set(["keys", "name", "modifiable", "calibration", "dim", "input_keys", "page"]))
        self.assertTrue(parameter_name(['condom','reg']) == 'Condom use proportion for regular sexual acts')

    def test_upload_data(self):
        import re
        import os
        import filecmp
        # create project
        response = self.client.post('/api/project/create/test', data = '{}')
        self.assertEqual(response.status_code, 200)
        # upload data
        example_excel_file_name = 'example.xlsx'
        file_path = helpers.safe_join(app.static_folder, example_excel_file_name)
        example_excel = open(file_path)
        headers = [('project', 'test')]
        response = self.client.post('api/project/update', headers=headers, data=dict(file=example_excel))
        example_excel.close()
        self.assertEqual(response.status_code, 200)
        # get data back and save the received file
        response = self.client.get('/api/project/workbook/test')
        content_disposition = response.headers.get('Content-Disposition')
        self.assertTrue(len(content_disposition)>0)
        file_name_info = re.search('filename=\s*(\S*)', content_disposition)
        self.assertTrue(len(file_name_info.groups())>0) 
        file_name = file_name_info.group(1)
        self.assertEqual(file_name,'test.xlsx')
        output_path = '/tmp/project_test.xlsx'
        if os.path.exists(output_path):os.remove(output_path)
        f = open(output_path, 'wb')
        f.write(response.data)
        f.close()
        # compare with source file
        result = filecmp.cmp(file_path, output_path)
        self.assertTrue(result)
        os.remove(output_path)

    def test_copy_project(self):
        # create project
        response = self.client.post('/api/project/create/test', data = '{}')
        self.assertEqual(response.status_code, 200)
        # upload data so that we can check the existence of data in the copied project
        example_excel_file_name = 'example.xlsx'
        file_path = helpers.safe_join(app.static_folder, example_excel_file_name)
        example_excel = open(file_path)
        headers = [('project', 'test')]
        response = self.client.post('api/project/update', headers=headers, data=dict(file=example_excel))
        example_excel.close()
        #get the info for the existing project
        response = self.client.get('/api/project/info', headers=headers)
        self.assertEqual(response.status_code, 200)
        old_info=json.loads(response.data)
        self.assertEqual(old_info['has_data'], True)
        response = self.client.post('/api/project/copy/test?to=test_copy', headers=headers)
        self.assertEqual(response.status_code, 200)
        #open the copy of the project
        response = self.client.get('/api/project/open/test_copy', headers=headers)
        self.assertEqual(response.status_code, 200)
        #get info for the copy of the project
        headers = [('project', 'test_copy')]
        response = self.client.get('/api/project/info', headers=headers)
        self.assertEqual(response.status_code, 200)
        new_info = json.loads(response.data)
        self.assertEqual(old_info['has_data'], True)
        #compare some elements
        self.assertEqual(old_info['populations'], new_info['populations'])
        self.assertEqual(old_info['programs'], new_info['programs'])


if __name__ == '__main__':
    unittest.main()
