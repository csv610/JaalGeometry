/****************************************************************************
** Meta object code from reading C++ file 'MeshConstraintsDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.5.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "MeshConstraintsDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MeshConstraintsDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.5.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JMeshConstraintsDialog_t {
    QByteArrayData data[15];
    char stringdata0[214];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JMeshConstraintsDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JMeshConstraintsDialog_t qt_meta_stringdata_JMeshConstraintsDialog = {
    {
        QT_MOC_LITERAL(0, 0, 22), // "JMeshConstraintsDialog"
        QT_MOC_LITERAL(1, 23, 17), // "setNewConstraints"
        QT_MOC_LITERAL(2, 41, 0), // ""
        QT_MOC_LITERAL(3, 42, 4), // "init"
        QT_MOC_LITERAL(4, 47, 22), // "addBoundaryConstraints"
        QT_MOC_LITERAL(5, 70, 17), // "invertConstraints"
        QT_MOC_LITERAL(6, 88, 16), // "setSelectionMode"
        QT_MOC_LITERAL(7, 105, 15), // "setSelectRegion"
        QT_MOC_LITERAL(8, 121, 17), // "openPolygonDialog"
        QT_MOC_LITERAL(9, 139, 10), // "getPolygon"
        QT_MOC_LITERAL(10, 150, 15), // "saveConstraints"
        QT_MOC_LITERAL(11, 166, 15), // "readConstraints"
        QT_MOC_LITERAL(12, 182, 10), // "lastDelete"
        QT_MOC_LITERAL(13, 193, 8), // "clearAll"
        QT_MOC_LITERAL(14, 202, 11) // "closeDialog"

    },
    "JMeshConstraintsDialog\0setNewConstraints\0"
    "\0init\0addBoundaryConstraints\0"
    "invertConstraints\0setSelectionMode\0"
    "setSelectRegion\0openPolygonDialog\0"
    "getPolygon\0saveConstraints\0readConstraints\0"
    "lastDelete\0clearAll\0closeDialog"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JMeshConstraintsDialog[] = {

// content:
    7,       // revision
    0,       // classname
    0,    0, // classinfo
    13,   14, // methods
    0,    0, // properties
    0,    0, // enums/sets
    0,    0, // constructors
    0,       // flags
    1,       // signalCount

// signals: name, argc, parameters, tag, flags
    1,    0,   79,    2, 0x06 /* Public */,

// slots: name, argc, parameters, tag, flags
    3,    0,   80,    2, 0x08 /* Private */,
    4,    0,   81,    2, 0x08 /* Private */,
    5,    0,   82,    2, 0x08 /* Private */,
    6,    0,   83,    2, 0x08 /* Private */,
    7,    0,   84,    2, 0x08 /* Private */,
    8,    0,   85,    2, 0x08 /* Private */,
    9,    0,   86,    2, 0x08 /* Private */,
    10,    0,   87,    2, 0x08 /* Private */,
    11,    0,   88,    2, 0x08 /* Private */,
    12,    0,   89,    2, 0x08 /* Private */,
    13,    0,   90,    2, 0x08 /* Private */,
    14,    0,   91,    2, 0x08 /* Private */,

// signals: parameters
    QMetaType::Void,

// slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

    0        // eod
};

void JMeshConstraintsDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JMeshConstraintsDialog *_t = static_cast<JMeshConstraintsDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0:
            _t->setNewConstraints();
            break;
        case 1:
            _t->init();
            break;
        case 2:
            _t->addBoundaryConstraints();
            break;
        case 3:
            _t->invertConstraints();
            break;
        case 4:
            _t->setSelectionMode();
            break;
        case 5:
            _t->setSelectRegion();
            break;
        case 6:
            _t->openPolygonDialog();
            break;
        case 7:
            _t->getPolygon();
            break;
        case 8:
            _t->saveConstraints();
            break;
        case 9:
            _t->readConstraints();
            break;
        case 10:
            _t->lastDelete();
            break;
        case 11:
            _t->clearAll();
            break;
        case 12:
            _t->closeDialog();
            break;
        default:
            ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (JMeshConstraintsDialog::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&JMeshConstraintsDialog::setNewConstraints)) {
                *result = 0;
            }
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JMeshConstraintsDialog::staticMetaObject = {
    {   &QDialog::staticMetaObject, qt_meta_stringdata_JMeshConstraintsDialog.data,
        qt_meta_data_JMeshConstraintsDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR
    }
};


const QMetaObject *JMeshConstraintsDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JMeshConstraintsDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JMeshConstraintsDialog.stringdata0))
        return static_cast<void*>(const_cast< JMeshConstraintsDialog*>(this));
    if (!strcmp(_clname, "Ui::MeshConstraintsDialog"))
        return static_cast< Ui::MeshConstraintsDialog*>(const_cast< JMeshConstraintsDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JMeshConstraintsDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 13)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 13;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 13)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 13;
    }
    return _id;
}

// SIGNAL 0
void JMeshConstraintsDialog::setNewConstraints()
{
    QMetaObject::activate(this, &staticMetaObject, 0, Q_NULLPTR);
}
QT_END_MOC_NAMESPACE
