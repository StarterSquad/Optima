define(['./../module', 'angular', 'underscore'], function (module, angular, _) {
  'use strict';

  module.controller('ProgramSetController', function ($scope, $http, programSetModalService,
    modalService, currentProject, projectApiService, $upload) {

    var openProject = currentProject.data;
    var defaults;
    var parameters;

    // Do not allow user to proceed if spreadsheet has not yet been uploaded for the project
    if (!openProject.has_data) {
      modalService.inform(
        function (){ },
        'Okay',
        'Please upload spreadsheet to proceed.',
        'Cannot proceed'
      );
      $scope.missingData = true;
      return;
    }

    // Get the list of saved programs from DB and set the first one as active
    $http.get('/api/project/' + openProject.id + '/progsets' )
      .success(function (response) {
        if(response.progsets) {
          console.log(JSON.stringify(response, null, 2));
          $scope.programSetList = response.progsets;
          if (response.progsets && response.progsets.length > 0) {
            $scope.activeProgramSet = response.progsets[0];
          }
        }
      });

    $scope.getCategories = function() {
      return _.uniq(_.pluck($scope.activeProgramSet.programs, "category"));
    }

    // Fetching default categories and programs for the open project
    projectApiService.getDefault(openProject.id)
      .success(function (response) {
        defaults = response;
      });

    // Get the list of default parameters for the project
    $http.get('/api/project/' + openProject.id + '/parameters' )
      .success(function (response) {
        parameters = response.parameters;
      });

    // This method is called by <select> to change current active program
    $scope.setActiveProgramSet = function(activeProgramSet) {
      $scope.activeProgramSet = activeProgramSet;
    };

    // Open pop-up to add new programSet
    $scope.addProgramSet = function () {
      var add = function (name) {
        var newProgramSet = {name:name, programs: angular.copy(defaults.programs)};
        $scope.programSetList[$scope.programSetList ? $scope.programSetList.length : 0] = newProgramSet;
        $scope.activeProgramSet = newProgramSet;
      };
      programSetModalService.openProgramSetModal(add, 'Add program set', $scope.programSetList, null, 'Add');
    };

    // Open pop-up to re-name programSet
    $scope.renameProgramSet = function () {
      if (!$scope.activeProgramSet) {
        modalService.informError([{message: 'No program set selected.'}]);
      } else {
        var rename = function (name) {
          $scope.activeProgramSet.name = name;
        };
        programSetModalService.openProgramSetModal(rename, 'Rename program set', $scope.programSetList, $scope.activeProgramSet.name, 'Update', true);
      }
    };

    // Download  programs data
    $scope.downloadProgramSet = function() {
      if(!$scope.activeProgramSet.id) {
        modalService.inform(
          function (){ },
          'Okay',
          'Please save the program set to proceed.',
          'Cannot proceed'
        );
      } else {
        $http.get('/api/project/' + openProject.id +  '/progsets' + '/' + $scope.activeProgramSet.id +'/data',
          {headers: {'Content-type': 'application/octet-stream'},
            responseType:'blob'})
          .success(function (response) {
            var blob = new Blob([response], { type: 'application/octet-stream' });
            saveAs(blob, ($scope.activeProgramSet.name + '.prj'));
          });
      }
    };

    // Upload programs data
    $scope.uploadProgramSet = function() {
      if(!$scope.activeProgramSet.id) {
        modalService.inform(
          function (){ },
          'Okay',
          'Please save the program set to proceed.',
          'Cannot proceed'
        );
      } else {
        angular
          .element('<input type=\'file\'>')
          .change(function (event) {
            $upload.upload({
              url: '/api/project/' + openProject.id + '/progsets' + '/' + $scope.activeProgramSet.id + '/data',
              file: event.target.files[0]
            }).success(function () {
              window.location.reload();
            });
          }).click();
      }
    };

    // Delete a programSet from $scope.programSetList and also from DB is it was saved.
    $scope.deleteProgramSet = function () {
      if (!$scope.activeProgramSet) {
        modalService.informError([{message: 'No program set selected.'}]);
      } else {
        var remove = function () {
          if ($scope.activeProgramSet.id) {
            $http.delete('/api/project/' + openProject.id +  '/progsets' + '/' + $scope.activeProgramSet.id).
              success(function() {
                deleteProgramSetFromPage();
              });
          } else {
            deleteProgramSetFromPage();
          }
        };
        modalService.confirm(
          function () {
            remove()
          }, function () {
          }, 'Yes, remove this program set', 'No',
          'Are you sure you want to permanently remove program set "' + $scope.activeProgramSet.name + '"?',
          'Delete program set'
        );
      }
    };

    var deleteProgramSetFromPage = function() {
      $scope.programSetList = _.filter($scope.programSetList, function (programSet) {
        return programSet.name !== $scope.activeProgramSet.name;
      });
      if($scope.programSetList && $scope.programSetList.length > 0) {
        $scope.activeProgramSet = $scope.programSetList[0];
      } else {
        $scope.activeProgramSet = undefined;
      }
    };

    // Copy a program-set
    $scope.copyProgramSet = function () {
      if (!$scope.activeProgramSet) {
        modalService.informError([{message: 'No program set selected.'}]);
      } else {
        var copy = function (name) {
          var copiedProgramSet = {name: name, programs: $scope.activeProgramSet.programs};
          $scope.programSetList[$scope.programSetList.length] = copiedProgramSet;
          $scope.activeProgramSet = copiedProgramSet;
        };
        programSetModalService.openProgramSetModal(copy, 'Copy program set', $scope.programSetList, $scope.activeProgramSet.name + ' copy', 'Copy');
      }
    };

    // Save a programSet to DB
    $scope.saveProgramSet = function() {
      var errorMessage;
      if (!$scope.activeProgramSet || !$scope.activeProgramSet.name) {
        errorMessage = 'Please create a new program set before trying to save it.';
      }
      if (errorMessage) {
        modalService.inform(
          function (){ },
          'Okay',
          errorMessage,
          'Cannot proceed'
        );
      } else {
        $http({
          url: '/api/project/' + openProject.id + '/progsets' + ($scope.activeProgramSet.id ? '/' + $scope.activeProgramSet.id : ''),
          method: ($scope.activeProgramSet.id ? 'PUT' : 'POST'),
          data: $scope.activeProgramSet
        }).success(function (response) {
          if(response.id) {
            $scope.activeProgramSet.id = response.id;
          }
          modalService.inform(
            function (){},
            'Okay',
            'Program set saved successfully',
            'Program set saved'
          );
        });
      }
    };

    // Opens a modal for editing an existing program.
    $scope.openEditProgramModal = function ($event, program) {
      var editProgram = angular.copy(program);
      editProgram.short_name = editProgram.short || editProgram.short_name;
      return programSetModalService.openProgramModal(editProgram, openProject, $scope.activeProgramSet.programs, parameters, $scope.getCategories()).result.then(
        function (newProgram) {
          $scope.activeProgramSet.programs[$scope.activeProgramSet.programs.indexOf(program)] = newProgram;
        }
      );
    };

    // Creates a new program and opens a modal for editing.
    $scope.openAddProgramModal = function ($event) {
      if ($event) {
        $event.preventDefault();
      }
      var program = {};

      return programSetModalService.openProgramModal(program, openProject, $scope.activeProgramSet.programs, parameters, $scope.getCategories()).result.then(
        function (newProgram) {
          $scope.activeProgramSet.programs.push(newProgram);
        }
      );
    };

    // Makes a copy of an existing program and opens a modal for editing.
    $scope.copyProgram = function ($event, existingProgram) {
      if ($event) {
        $event.preventDefault();
      }
      var program = angular.copy(existingProgram);
      program.name = program.name + ' copy';
      program.short_name = (program.short || program.short_name ) + ' copy';
      program.short = undefined;

      return programSetModalService.openProgramModal(program, openProject, $scope.activeProgramSet.programs, parameters, $scope.getCategories()).result.then(
        function (newProgram) {
          $scope.activeProgramSet.programs.push(newProgram);
        }
      );
    };

  });
});
