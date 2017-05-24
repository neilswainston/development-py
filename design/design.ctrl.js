designApp.controller("designCtrl", function($scope) {

    $scope.models = {
        selected: null,
        templates: [
        	{type: "design", id: 2, "features": []},
            {type: "feature", id: 1}
        ],
        designs: [
        	{
        		type: "design",
        		id: 1, "features": []
        	}
        ]
    };

    $scope.$watch('models.designs', function(model) {
        $scope.modelAsJson = angular.toJson(model, true);
    }, true);
});