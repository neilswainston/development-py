designApp.directive("designPanel", function() {
    return {
    	scope: {
    		"addDesign": "&",
    		"removeDesign": "&",
    		"models": "="
    	},
        templateUrl: "design.html"
    };
});