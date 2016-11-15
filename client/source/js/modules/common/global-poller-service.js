/**
 * poller is a factory for tracking jobs on the server
 */

define(['angular' ], function (angular) {

  'use strict';

  return angular
    .module('app.common.global-poller-service', [])
    .factory('globalPoller', ['$http', '$timeout', function($http, $timeout) {

      var polls = {};

      function getPoll(id) {
        if (!(id in polls)) {
          console.log('Creating polling slot for', id);
          polls[id] = {isRunning: false, id: id};
        }
        return polls[id];
      }

      function startPoll(id, url, callback) {
        var poll = getPoll(id);
        poll.url = url;
        poll.callback = callback;

        if (!poll.isRunning) {
          console.log('Launch polling for', poll.id);
          poll.isRunning = true;

          function pollWithTimeout() {
            var poll = getPoll(id);
            $http
              .get(poll.url)
              .success(function(response) {
                console.log(response);
                if (response.status === 'started') {
                  poll.timer = $timeout(pollWithTimeout, 1000);
                } else {
                  stopPoll(id);
                }
                poll.callback(response);
              })
              .error(function(response) {
                stopPoll(id);
                poll.callback(response);
              });
          }

          pollWithTimeout();
        }
      }

      function stopPoll(id) {
        var poll = getPoll(id);
        if (poll.isRunning) {
          console.log('Stop polling for', poll.id);
          poll.isRunning = false;
          $timeout.cancel(poll.timer);
        };
      }

      function stopPolls() {
        _.each(polls, function(poll) {
          if (poll.isRunning) {
            console.log('Stop polling for', poll.id);
            poll.isRunning = false;
            $timeout.cancel(poll.timer);
          }
        });
      }

      function killJob(projectId, workType) {
        $http
          .delete('/api/task/' + projectId + '/type/' + workType);
      }

      return {
        startPoll: startPoll,
        stopPolls: stopPolls,
        killJob: killJob
      };

    }]);


});