import React, {useContext} from 'react';
import {useState} from "react";
import {useNavigate} from "react-router-dom";

export default function ToolDescription({oid, URL, imagePath, toolTitle, textDescription, isFavoritable=true}) {

    let favoritesList = localStorage.getItem('favorites')
    const [isFavorited, setIsFavorited] = useState(
        favoritesList !== null ? favoritesList.indexOf(oid) >= 0 : false
    );



    const onClickHeart = (e) => {
        e.stopPropagation();
        setIsFavorited(!isFavorited);

        //adds or removes favorites from localStorage
        let newFavoriteArray = JSON.parse(localStorage.getItem('favorites'));
        if (isFavorited) { //value of isFavorited is from before toggle, since it only gets updated after render
            let index = newFavoriteArray.indexOf(oid);
            newFavoriteArray.splice(index, 1);
        }
        else {
            newFavoriteArray.push(oid)
        }
        localStorage.setItem('favorites', JSON.stringify(newFavoriteArray));

        updateToolPopularity();
    }



    //updates the tool's popularity on the server
    function updateToolPopularity() {
        fetch("http://localhost:5000/update-tool-popularities", {
            method: "POST",
            headers: {'Content-Type': 'application/json'},
            body: JSON.stringify({"oid": oid,"isFavorited": isFavorited})})
            .then(response => {
                if (response.ok) {
                    console.log('successfully updated tool popularity');
                } else {
                    console.log('failed to update tool popularity');
                }
            })
            .catch(error => {
                console.error('Error:', error);
            });
    }



    const heartIconLink = isFavorited ? '/images/filled_in_heart.png' : '/images/empty_heart.png'

    const navigate = useNavigate();
    const navigateToToolPage = () => {
        navigate(`/tool-page/${URL}`);
    }

    return (
      <div className={'container tool-description'} onClick={navigateToToolPage}>
          {isFavoritable ? <img className={'heart-icon'} src={heartIconLink} alt={'heart'} onClick={onClickHeart}/> : ''}
          <img className={'tool-reference-image'} src={imagePath} alt={'description'}/>
          <h3 className={'tool-title'}>{toolTitle}</h3>
          <p className={'tool-text-description'}>{textDescription}</p>
      </div>
    );
}