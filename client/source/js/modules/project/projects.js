define(['angular', 'ui.router'], function (angular) {

  'use strict';

  var module = angular.module('app.project', ['ui.router']);


  module.config(function ($stateProvider) {
    $stateProvider
      .state('home', {
        url: '/',
        templateUrl: 'js/modules/project/manage-projects.html?cacheBust=xxx',
        controller: 'ProjectOpenController'
      })
  });


  module.controller(
    'ProjectOpenController',
    function($scope, rpcService, modalService, userManager, projectService, $state, $upload, $modal, toastr) {

      function initialize() {
        $scope.sortType = 'name'; // set the default sort type
        $scope.sortReverse = false;  // set the default sort order
        $scope.projectService = projectService;
      }

      function getProjectNames() {
        return _.pluck(projectService.projects, 'name');
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
        _.forEach(projectService.projects, function(project) {
          project.selected = $scope.allSelected;
        });
      };

      $scope.deleteSelected = function() {
        const projectIds =
          _.filter(
            projectService.projects,
            function(project) { return project.selected; })
          .map(
            function(project) { return project.id; });
        projectService.deleteProjects(projectIds)
      };

      $scope.downloadSelected = function() {
        const projectsIds =
          _.filter(
            projectService.projects,
            function(project) { return project.selected; })
          .map(
            function(project) { return project.id; });
        projectService
          .downloadSelectedProjects(projectsIds)
          .success(function (response) {
            saveAs(new Blob([response], { type: "application/octet-stream", responseType: 'arraybuffer' }), 'portfolio.zip');
          });
      };

      $scope.open = function (name, id) {
        projectService.setActiveProjectId(id);
        projectService.getActiveProject()
          .then(function() {
            toastr.success('Project "' + name + '" loaded');
          });
      };

      $scope.copy = function(name, projectId) {
        projectService
          .copyProject(
            projectId,
            rpcService.getUniqueName(name, getProjectNames()))
          .then(function() {
            toastr.success('Copied project');
            $state.reload();
          });
      };

      $scope.copyOptimaLiteProject = function(project) {
        var name =
        projectService
          .copyProject(
            project.id,
            rpcService.getUniqueName(project.name, getProjectNames()))
          .then(function() {
            toastr.success('Project "'+project.name+'" loaded from database. Please proceed directly to analysis (scenarios and/or optimizations)');
            $state.reload();
          });
      };

      $scope.downloadSpreadsheet = function (name, id) {
        rpcService.rpcDownload(
          'download_data_spreadsheet', [id], {'is_blank': false})
        .then(function (response) {
          toastr.success('Spreadsheet downloaded');
        });
      };

      $scope.uploadProject = function() {
        projectService
          .uploadProject()
          .then(function() {
            toastr.success('Project uploaded');
          });
      };

      $scope.uploadProjectFromSpreadsheet = function() {
        projectService
          .uploadProjectFromSpreadsheet()
          .then(function() {
            toastr.success('Project uploaded from spreadsheet');
          });
      };

      $scope.uploadSpreadsheet = function(projectName, projectId) {
        rpcService
          .rpcUpload(
            'update_project_from_uploaded_spreadsheet', [projectId])
          .then(function(response) {
            toastr.success('Uploaded spreadsheet for project');
          });
      };

      $scope.editProjectName = function(project) {
        modalService.rename(
          function(name) {
            project.name = name;
            projectService
              .renameProject(project.id, project)
              .then(function () {
                toastr.success('Renamed project');
                $state.reload();
              });
          },
          'Edit project name',
          "Enter project name",
          project.name,
          "Name already exists",
          _.without(getProjectNames(), project.name));
      };

      $scope.downloadProject = function (name, id) {
        rpcService
          .rpcDownload(
            'download_project', [id])
          .then(function() {
            toastr.success('Project downloaded');
          });
      };

      $scope.downloadPrjWithResults = function (name, id) {
        rpcService
          .rpcDownload(
            'download_project_with_result', [id])
          .then(function() {
            toastr.success('Project downloaded');
          });
      };

      $scope.openOptimaLiteProjectList = function() {
        modalService.optimaLiteProjectList();
      };

      initialize();
  });

  return module;

});

