designApp.directive("designPanel", function() {
    return {
    	scope: {
    		templates: "=",
    		query: "&",
    		selected: "&",
    		setSelected: "&",
    		addDesign: "&",
    		removeDesign: "&",
    	},
        templateUrl: "design.html"
    };
});