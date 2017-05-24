designApp.controller("designCtrl", function($scope) {

    $scope.models = {
        selected: null,
        templates: [
        	{type: "design", typeName: "Design", id: 2, "features": []},
            {type: "feature", typeName: "Feature", id: 1}
        ],
        designs: [
        	{
        		type: "design",
        		typeName: "Design",
        		id: 1, "features": []
        	}
        ]
    };

    $scope.$watch('models.designs', function(model) {
        $scope.modelAsJson = angular.toJson(model, true);
    }, true);
});