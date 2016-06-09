define(['./module', 'angular', 'underscore'], function (module, angular, _) {
  'use strict';

  module.controller('ModelCalibrationController', function ($scope, $http, info, modalService, $upload, $modal, $timeout) {

    var activeProjectInfo = info.data;
    var defaultParameters;
    $scope.parsets = [];
    $scope.activeParset = undefined;
    $scope.state = {maxtime: '10'};

    // Check if current active project has spreadsheet uploaded for it.
    if (!activeProjectInfo.has_data) {
      modalService.inform(
        function (){ },
        'Okay',
        'Please upload spreadsheet to proceed.',
        'Cannot proceed'
      );
      $scope.missingData = true;
      return;
    }

    // Fetching list of parsets for open project
    $http.get('/api/project/' + activeProjectInfo.id + '/parsets').
      success(function (response) {
        var parsets = response.parsets;
        if(parsets) {
          $scope.parsets = parsets;
          $scope.activeParset = parsets[0];
          $scope.getGraphs();
        }
      });

    // Fetching graphs for active parset
    $scope.getGraphs = function() {
      var data = {};
      if($scope.selectors) {
        var selectors = _.filter($scope.selectors, function(selector) {
          return selector.checked;
        }).map(function(selector) {
          return selector.key;
        });
        if(selectors && selectors.length > 0) {
          data.which = selectors;
        }
      }
      $http.get('/api/project/' + activeProjectInfo.id + '/parsets/' + $scope.activeParset.id + '/calibration', {params: data})
        .success(function (response) {
          setCalibrationData(response.calibration);
          // console.log(JSON.stringify(response.calibration.parameters, null, 2))
        });
    };

    // Sending parameters to re-process graphs for active parset
    $scope.processGraphs = function(shouldSave) {
      var data = {};
      if($scope.parameters) {
        data.parameters = $scope.parameters;
      }
      if($scope.selectors) {
        var selectors = _.filter($scope.selectors, function(selector) {
          return selector.checked;
        }).map(function(selector) {
          return selector.key;
        });
        if(selectors && selectors.length > 0) {
          data.which = selectors;
        }
      }
      var url = '/api/project/' + activeProjectInfo.id + '/parsets/' + $scope.activeParset.id + '/calibration';
      if (shouldSave) {
        url = url + '?doSave=true';
      }
      $http.put(url, data)
      .success(function (response) {
        setCalibrationData(response.calibration);
      });
    };

    // Set calibration data in scope
    var setCalibrationData = function(calibration) {
      console.log(calibration);
      $scope.graphs = calibration.graphs;
      // $scope.calibrationChart = calibration.graphs;
      // $scope.selectors = calibration.selectors;
      defaultParameters = calibration.parameters;
      $scope.parameters = angular.copy(calibration.parameters);
    };

    // Reset changes made in parameters
    $scope.resetParameters = function() {
      $scope.parameters = angular.copy(defaultParameters);
    };

    // Add parameter set
    $scope.addParameterSet = function() {
      var add = function (name) {
        $http.post(
            '/api/project/' + activeProjectInfo.id + '/parsets',
            { name: name})
        .success(function(response) {
          $scope.parsets = response;
          $scope.activeParset = response[response.length - 1];
        });
      };
      openParameterSetModal(add, 'Add parameter set', $scope.parsets, null, 'Add');
    };

    // Copy parameter set
    $scope.copyParameterSet = function() {
      if (!$scope.activeParset) {
        modalService.informError([{message: 'No parameter set selected.'}]);
      } else {
        var rename = function (name) {
          $http.post('/api/project/' + activeProjectInfo.id + '/parsets', {
            name: name,
            parset_id: $scope.activeParset.id
          }).success(function (response) {
            $scope.parsets = response;
            $scope.activeParset = response[response.length - 1];
          });
        };
        openParameterSetModal(rename, 'Copy parameter set', $scope.parsets, $scope.activeParset.name + ' copy', 'Copy');
      }
    };

    // Rename parameter set
    $scope.renameParameterSet = function() {
      if (!$scope.activeParset) {
        modalService.informError([{message: 'No parameter set selected.'}]);
      } else {
        var rename = function (name) {
          $http.put('/api/project/' + activeProjectInfo.id + '/parsets/' + $scope.activeParset.id, {
            name: name
          }).success(function () {
            $scope.activeParset.name = name;
          });
        };
        openParameterSetModal(rename, 'Rename parameter set', $scope.parsets, $scope.activeParset.name, 'Rename', true);
      }
    };

    // Delete parameter set
    $scope.deleteParameterSet = function() {
      if (!$scope.activeParset) {
        modalService.informError([{message: 'No parameter set selected.'}]);
      } else {
        var remove = function () {
          $http.delete('/api/project/' + activeProjectInfo.id + '/parsets/' + $scope.activeParset.id)
            .success(function() {
              $scope.parsets = _.filter($scope.parsets, function (parset) {
                return parset.id !== $scope.activeParset.id;
              });
              if($scope.parsets.length > 0) {
                $scope.activeParset = $scope.parsets[0];
              }
            });
        };
        // This has been temporarily commented out: https://trello.com/c/omuvJSYD/853-reuploading-spreadsheets-fails-in-several-ways
        // if ($scope.activeParset.name === "default") {
        if ( false ) {
          modalService.informError([{message: 'Deleting the default parameter set is not permitted.'}]);
        } else {
          modalService.confirm(
            function () {
              remove()
            }, function () {
            }, 'Yes, remove this parameter set', 'No',
            'Are you sure you want to permanently remove parameter set "' + $scope.activeParset.name + '"?',
            'Delete parameter set'
          );
        }
      }
    };

    // Download  parameter-set data
    $scope.downloadParameterSet = function() {
      $http.get('/api/project/' + activeProjectInfo.id +  '/parsets' + '/' + $scope.activeParset.id +'/data',
        {headers: {'Content-type': 'application/octet-stream'},
          responseType:'blob'})
        .success(function (response) {
          var blob = new Blob([response], { type: 'application/octet-stream' });
          saveAs(blob, ($scope.activeParset.name + '.par'));
        });
    };

    // Upload parameter-set data
    $scope.uploadParameterSet = function() {
      angular
        .element('<input type=\'file\'>')
        .change(function(event){
          $upload.upload({
            url: '/api/project/' + activeProjectInfo.id +  '/parsets/' + $scope.activeParset.id + '/data',
            file: event.target.files[0]
          }).success(function () {
            window.location.reload();
          });
        }).click();
    };

    // Opens modal to add / rename / copy parameter set
    var openParameterSetModal = function (callback, title, parameterSetList, parameterSetName, operation, isRename) {
      var onModalKeyDown = function (event) {
        if(event.keyCode == 27) { return modalInstance.dismiss('ESC'); }
      };
      var modalInstance = $modal.open({
        templateUrl: 'js/modules/model/parameter-set-modal.html',
        controller: ['$scope', '$document', function ($scope, $document) {
          $scope.title = title;
          $scope.name = parameterSetName;
          $scope.operation = operation;
          $scope.updateParameterSet = function () {
            $scope.newParameterSetName = $scope.name;
            callback($scope.name);
            modalInstance.close();
          };
          $scope.isUniqueName = function (parameterSetForm) {
            var exists = _(parameterSetList).some(function(item) {
                return item.name == $scope.name;
              }) && $scope.name !== parameterSetName && $scope.name !== $scope.newParameterSetName;
            if(isRename) {
              parameterSetForm.parameterSetName.$setValidity("parameterSetUpdated", $scope.name !== parameterSetName);
            }
            parameterSetForm.parameterSetName.$setValidity("parameterSetExists", !exists);
            return exists;
          };
          $document.on('keydown', onModalKeyDown); // observe
          $scope.$on('$destroy', function (){ $document.off('keydown', onModalKeyDown); });  // unobserve
        }]
      });
      return modalInstance;
    }

    $scope.startAutoCalibration = function() {
      var data = {};
      if ($scope.state.maxtime) {
        data.maxtime = Number($scope.state.maxtime);
      }
      $http.post('/api/project/' + activeProjectInfo.id +  '/parsets' + '/' + $scope.activeParset.id +'/automatic_calibration',
        data
      )
        .success(function(response) {
          if(response.status === 'started') {
            $scope.statusMessage = 'Automatic calibration started.';
            $scope.secondsRun = 0;
            $scope.setMaxtime = data.maxtime;
            pollAutoCalibration();
          } else if(response.status === 'running') {
            $scope.statusMessage = 'Automatic calibration already running.'
          }
        })
    };

    var pollAutoCalibration = function() {
      $http.get('/api/project/' + activeProjectInfo.id +  '/parsets/' + $scope.activeParset.id +'/automatic_calibration')
        .success(function(response) {
          if(response.status === 'completed') {
            getAutoCalibratedGraphs();
            $scope.statusMessage = 'Automatic calibration successfully completed.';
            $timeout.cancel($scope.pollTimer);
          } else if(response.status === 'started'){
            $scope.pollTimer = $timeout(pollAutoCalibration, 1000);
            $scope.statusMessage = "Running: " + $scope.secondsRun + " s.";
            $scope.secondsRun += 1;
          }
        });
    };

    var getAutoCalibratedGraphs = function() {
      $http.get(
        '/api/project/' + activeProjectInfo.id
        + '/parsets' + '/' + $scope.activeParset.id
        + '/calibration?autofit=true')
      .success(function(response) {
        setCalibrationData(response.calibration);
      });
    };

    $scope.resetAutoCalibration = function() {
      $scope.processGraphs();
    };

    $scope.saveAutoCalibration = function() {
      $http.post('/api/project/' + activeProjectInfo.id + '/parsets/' + $scope.activeParset.id + '/calibration', {
        parameters: $scope.parameters,
        result_id: $scope.result_id
      }).success(function(response) {
        $scope.statusMessage = 'Result saved';
      });
    };

  });
});
