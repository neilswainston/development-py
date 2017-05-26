designApp.factory("DesignService", [function() {
	var obj = {};
	
	obj.query = {
		"app": "PartsGenie",
		"designs": [], 
		"filters": {
			"max_repeats": 5
		},
	};
	
	obj.selected = null;
	
	obj.setSelected = function(selected) {
		obj.selected = selected;
	}

	return obj;
}]);