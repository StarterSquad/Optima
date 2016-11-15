define([
  'angular',
  'ng-resource'
], function (angular) {

  angular.module('app.resources.model', [
    'ngResource'
  ])
    .service('Model', function ($resource) {
      return $resource('/api/model/:path/:suffix/:postsuffix',
        { path: '@path' },
        {
          runManualCalibration: {
            method: 'POST',
            isArray: false,
            params: {
              path: 'calibrate',
              suffix: 'manual'
            }
          },
          getKeyDataMeta: {
            method: 'GET',
            isArray: false,
            params: {
              path: 'data',
              suffix: 'data',
              postsuffix: 'meta'
            }
          },
          getCalibrateParameters: {
            method: 'GET',
            isArray: false,
            params: {
              path: 'calibrate',
              suffix: 'parameters'
            }
          },
          getPrograms: {
            method: 'GET',
            isArray: false,
            params: {
              path: 'data',
              suffix: 'programs'
            }
          }
        }
      );
    });
});
