define(['./module', 'angular', 'underscore'], function (module, angular, _) {
  'use strict';

  module.controller('ProjectOpenController',
    function ($scope, $http, activeProject, projects, modalService,
              fileUpload, UserManager, projectApiService, $state, toastr) {

      function initialize() {
        $scope.sortType = 'name'; // set the default sort type
        $scope.sortReverse = false;  // set the default sort order
        $scope.activeProjectId = activeProject.getProjectIdForCurrentUser();
        loadProjects(projects.data.projects);
      }

      function loadProjects(projects) {
        $scope.projects = _.map(projects, function(project) {
          project.creationTime = Date.parse(project.creationTime);
          project.updatedTime = Date.parse(project.updatedTime);
          project.dataUploadTime = Date.parse(project.dataUploadTime);
          return project;
        });
        console.log('projects', $scope.projects);
      }

      $scope.filterByName = function(project) {
        if ($scope.searchTerm) {
          return project.name.toLowerCase().indexOf($scope.searchTerm.toLowerCase()) !== -1;
        }
        return true;
      };

      $scope.updateSorting = function(sortType) {
        if ($scope.sortType === sortType) {
          $scope.sortReverse = !$scope.sortReverse;
        } else {
          $scope.sortType = sortType;
        }
      };

      $scope.selectAll = function() {
        _.forEach($scope.projects, function(project) {
          project.selected = $scope.allSelected;
        });
      };

      $scope.deleteSelected = function() {
        const selectedProjectIds = _.filter($scope.projects, function(project) {
          return project.selected;
        }).map(function(project) {
          return project.id;
        });
        projectApiService.deleteSelectedProjects(selectedProjectIds)
          .success(function () {
            $scope.projects = _.filter($scope.projects, function(project) {
              return !project.selected;
            });
            _.each(selectedProjectIds, function(projectId) {
              activeProject.ifActiveResetFor(projectId, UserManager.data);
            });
          });
      };

      $scope.downloadSelected = function() {
        const selectedProjectsIds = _.filter($scope.projects, function(project) {
          return project.selected;
        }).map(function(project) {
          return project.id;
        });
        projectApiService.downloadSelectedProjects(selectedProjectsIds)
          .success(function (response) {
            saveAs(new Blob([response], { type: "application/octet-stream", responseType: 'arraybuffer' }), 'portfolio.zip');
          });
      };

      $scope.open = function (name, id) {
        activeProject.setActiveProjectFor(name, id, UserManager.data);
        $state.go('home');
      };

      function getUniqueName(name, otherNames) {
        var i = 0;
        var uniqueName = name;
        while (_.indexOf(otherNames, uniqueName) >= 0) {
          i += 1;
          uniqueName = name + ' (' + i + ')';
        }
        return uniqueName;
      }

      $scope.copy = function(name, id) {
        var otherNames = _.pluck($scope.projects, 'name');
        var newName = getUniqueName(name + '(Copy)', otherNames);
        projectApiService.copyProject(id, newName).success(function (response) {
          projectApiService.getProjectList()
            .success(function(response) {
              toastr.success('Copied project ' + newName);
              loadProjects(response.projects);
            });
        });
      };

      /**
       * Opens to edit an existing project using name and id in /project/create screen.
       */
      $scope.edit = function (name, id) {
        activeProject.setActiveProjectFor(name, id, UserManager.data);
        $state.go('project.edit');
      };

      /**
       * Regenerates workbook for the given project.
       */
      $scope.workbook = function (name, id) {
        // read that this is the universal method which should work everywhere in
        // http://stackoverflow.com/questions/24080018/download-file-from-a-webapi-method-using-angularjs
        window.open(projectApiService.getSpreadsheetUrl(id), '_blank', '');
      };

      /**
       * Gets the data for the given project `name` as <name>.json  file.
       */
      $scope.getData = function (name, id) {
        projectApiService.getProjectData(id)
          .success(function (response, status, headers, config) {
            var blob = new Blob([response], { type: 'application/octet-stream' });
            saveAs(blob, (name + '.prj'));
          });
      };

      /**
       * Upload data spreadsheet for a project.
       */
      $scope.setData = function (name, id, file) {
        var message = 'Warning: This will overwrite ALL data in the project ' + name + '. Are you sure you wish to continue?';
        modalService.confirm(
          function (){ fileUpload.uploadDataSpreadsheet($scope, file, projectApiService.getDataUploadUrl(id), false); },
          function (){},
          'Yes, overwrite data',
          'No',
          message,
          'Upload data'
        );
      };

      /**
       * Upload project data.
       */
      $scope.preSetData = function(name, id) {
        angular
          .element('<input type=\'file\'>')
          .change(function(event){
          $scope.setData(name, id, event.target.files[0]);
        }).click();
      };

      /**
       * Removes the project.
       */
      var removeProject = function (name, id, index) {
        projectApiService.deleteProject(id).success(function (response) {
          $scope.projects = _($scope.projects).filter(function (item) {
            return item.id != id;
          });
          activeProject.ifActiveResetFor(id, UserManager.data);
        });
      };

      /**
       * Opens a modal window to ask the user for confirmation to remove the project and
       * removes the project if the user confirms.
       * Closes it without further action otherwise.
       */
      $scope.remove = function ($event, name, id, index) {
        if ($event) { $event.preventDefault(); }
        var message = 'Are you sure you want to permanently remove project "' + name + '"?';
        modalService.confirm(
          function (){ removeProject(name, id, index); },
          function (){  },
          'Yes, remove this project',
          'No',
          message,
          'Remove project'
        );
      };

      initialize();
  });

});
